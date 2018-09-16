/* Copyright (C) 2013-2016, The Regents of The University of Michigan.
All rights reserved.

This software was developed in the APRIL Robotics Lab under the
direction of Edwin Olson, ebolson@umich.edu. This software may be
available under alternative licensing terms; contact the address above.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the Regents of The University of Michigan.
*/

#include "beacon_calibration.h"

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#define SPEED_OF_LIGHT 299.792458               // meters per microsecond

struct calib_data {
    int n_locs;
    int n_delays;
    int n_tprops;
    bool disable_z; // makes 3d behave like 2d by effectively keeping z at 0
    // delay index to use to define the z-plane. The intention here is for the id to correspond
    // to the robot beacon which performs the "calibration dance," moving over the z plane
    // and collecting data at 4-6 points, which we can then all set to be z=0.
    int z_plane_delay_idx;
    // It may be useful in the 2d case to still have the z-plane defined by z_plane_delay_idx
    // at a different z level than the rest of the beacons. Say because they are on the floor
    // and the robot is at a different level. It is still 2d calibration because z coords are
    // not modified at all.
    double z_plane_2d_offset;
    zarray_t **tprop_peaks;
    measurement_density_t *densities;
    double *initial_tprops;
    int *as; // location index "a" for each tprop, size n_tprops
    int *bs; // location index "b" for each tprop, size n_tprops
    int *ds; // antenna delay index for each location index, size n_locs
    double *coords; // size 2 * n_locs
    double *delays;
    bool opt_coords;
    bool opt_delays;
    double expected_delay;
    double delay_sigma;

    // selects out the peak indices to be used by dist_for in multilateration,
    // so that different combinations may be attempted
    int *multilateration_peak_choices;

    // enables diagnostics to evaluate which peaks the max mixture chooses
    bool count_peak_choices;
    int chosen_peaks[16];
};

typedef struct calib_data calib_data_t;

// based on https://www.ebayinc.com/stories/blogs/tech/fast-approximate-logarithms-part-iii-the-formulas/
double fastlog2_base(double x, double a, double b, double c)    // compute log2(x) by reducing x to [0.75, 1.5)
{
    /** MODIFY THIS SECTION **/
    // (x-1)*(a*(x-1) + b)/((x-1) + c) (line 8 of table 2)
    #define FN (fexp + signif*(a*signif + b)/(signif + c))

    union { double f; uint64_t i; } ux1, ux2;
    /*
     * Assume IEEE representation, which is sgn(1):exp(11):frac(52)
     * representing (1+frac)*2^(exp-1023).  Call 1+frac the significand
     */

    // get exponent
    ux1.f = x;
    int iexp = (ux1.i >> 52L) & 0x07ff;

    int greater = ux1.i & (1L << 51);  // true if signifand > 1.5
    if (greater) {
        // significand >= 1.5 so need to divide by 2.  Accomplish this by
        // stuffing exp = 1022 which corresponds to an exponent of -1
        ux2.i = (ux1.i & 0x0fffffffffffffL) | 0x3fe0000000000000L;
        // 1022 instead of 1023 compensates for division by 2
        double fexp = iexp - 1022;
        double signif = ux2.f - 1.0;
        return FN;
    } else {
        // get significand by stuffing exp = 1023 which corresponds to an exponent of 0
        ux2.i = (ux1.i & 0x0fffffffffffffL) | 0x3ff0000000000000L;
        double fexp = iexp - 1023;
        double signif = ux2.f - 1.0;
        return FN;
    }
    #undef FN
}

double fastlog2(double x)
{
    // minimizes max relative error: 0.326595350590, 2.198627976298, 1.523722925329
    // minimizes mean square relative error : 0.271489684597, 2.380464602029, 1.650898108180
    return fastlog2_base(x, 0.271489684597, 2.380464602029, 1.650898108180);
}

double fastlog(double x)
{
    return fastlog2(x) * 0.693147180559945309417232; //  * log(2)
}

// based on http://www.machinedlearnings.com/2011/06/fast-approximate-logarithm-exponential.html
double fastpow2_base(double p, double a, double b, double c, double d)
{
    // prevent overflow onto sign later
    if (p < -1022) {
        return 0;
    }

    union { double f; uint64_t i; } vp = { p };
    int sign = (vp.i >> 63);
    double z = p - (int)p + sign;
    union { uint64_t i; double f; } v = { (uint64_t)((1L << 52) * (p + a + b / (c - z) - d * z)) };
    return v.f;
}

double fastpow2(double p)
{
    return fastpow2_base(p, 1017.254358117020, 27.866401707383, 4.849981542817, 1.492398782372);
}

double fastexp(double p)
{
    return fastpow2(1.4426950408889634073599 * p); // log2(e) * p
}

// in meters
static double bias_by_signal_strength[17] = {
    -0.110, -0.105, -0.100, -0.093, -0.082, -0.069, -0.051, -0.027, 0.000,
    0.021, 0.035, 0.042, 0.049, 0.062, 0.071, 0.076, 0.081
};

static double bias_quadratic_as[17];
static double bias_quadratic_bs[17];
static double bias_quadratic_cs[17];

void initialize_range_bias_quadratic_fit()
{
    int n = 17;

    // the dx has been normalized to 1 and min_x to 0
    double inv_4_dx2 = 1.0 / 4.0;
    double inv_2_dx = 1.0 / 2.0;

    for (int i = 1; i < n - 1; i++) {
        double x2 = i;
        double y1 = bias_by_signal_strength[i - 1];
        double y2 = bias_by_signal_strength[i];
        double y3 = bias_by_signal_strength[i + 1];
        double a = (y3 - 2 * y2 + y1) * inv_4_dx2;
        double b = (y3 - y1) * inv_2_dx - a * 2.0 * x2;
        double c = y2 - a * x2 * x2 - b * x2;

        bias_quadratic_as[i] = a;
        bias_quadratic_bs[i] = b;
        bias_quadratic_cs[i] = c;
    }

    bias_quadratic_as[0] = bias_quadratic_as[1];
    bias_quadratic_bs[0] = bias_quadratic_bs[1];
    bias_quadratic_cs[0] = bias_quadratic_cs[1];

    bias_quadratic_as[n - 1] = bias_quadratic_as[n - 2];
    bias_quadratic_bs[n - 1] = bias_quadratic_bs[n - 2];
    bias_quadratic_cs[n - 1] = bias_quadratic_cs[n - 2];
}

// For Channel 1, 500MHz bandwidth, 64MHz PRF
// Application note, aps011_sources_of_error_in_twr
double range_bias_at(double m)
{
    if (!isfinite(m)) {
        fprintf(stderr, "Warning: NaN/Inf range value\n");
        exit(1);
        return 0;
    }
    if (m <= 0) {
        return 0;
    }
    double range_log_part = 10.6425824807107 + fastlog2(m) * 0.3010299995664; // log10((double)(4 * M_PI * 3494.4e6 * m));

    double signal_strength = -14.3 + 169.536 - 20 * range_log_part;
    double index_loc = (-61 - signal_strength) * 0.5;
    if (index_loc <= 0.0) {
        return bias_by_signal_strength[0];
    }

    int index_low = (int)index_loc;
    if (index_low >= 16) {
        return bias_by_signal_strength[16];
    }
    //
    // int rank = 12;
    // double coefs[] = { -7.26287179e-12,  5.61954307e-10, -1.76235568e-08,  2.76547983e-07, -1.95404578e-06, -2.39883414e-06,  1.52542833e-04, -1.19491591e-03, 4.49414628e-03, -8.46626765e-03,  7.39846074e-03,  2.64537844e-03, -1.10002282e-01 };
    // double bias = coefs[rank];
    // double factor = index_loc;
    // for (int i = 0; i < rank; i++) {
    //     bias += factor * coefs[rank - 1 - i];
    //     factor *= index_loc;
    // }

    if ("quadratic") {
        int best_i = (int)(index_loc + 0.5);
        double a = bias_quadratic_as[best_i];
        double b = bias_quadratic_bs[best_i];
        double c = bias_quadratic_cs[best_i];
        double bias = a * index_loc * index_loc + b * index_loc + c;
        return bias;
    } else {
        double fraction_high = index_loc - index_low;
        double bias = bias_by_signal_strength[index_low] * (1.0 - fraction_high) +
                      bias_by_signal_strength[index_low + 1] * fraction_high;
        return bias;
    }
}

double dist_from_tof(double flight_time)
{
    double dist = flight_time * SPEED_OF_LIGHT;
    return dist - range_bias_at(dist);
}

double dist_between_simple(double *coords, int a, int b)
{
    double dx = coords[2 * a] - coords[2 * b];
    double dy = coords[2 * a + 1] - coords[2 * b + 1];
    return sqrt(dx * dx + dy * dy);
}

double dist_between(calib_data_t *data, double *coords, int a, int b)
{
    double dx = coords[2 * a] - coords[2 * b];
    double dy = coords[2 * a + 1] - coords[2 * b + 1];
    double z1 = (data->ds[a] == data->z_plane_delay_idx) ? data->z_plane_2d_offset : 0;
    double z2 = (data->ds[b] == data->z_plane_delay_idx) ? data->z_plane_2d_offset : 0;
    double dz = z1 - z2;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    return dist;
}

double dist_for(calib_data_t *data, int a, int b) {
    if (a > b) {
        int tmp = a;
        a = b;
        b = tmp;
    }
    int idx = -1;
    for (int i = 0; i < data->n_tprops; i++) {
        if ((data->as[i] == a && data->bs[i] == b) ||
            (data->as[i] == b && data->bs[i] == a)) {
            idx = i;
            break;
        }
    }
    if (idx == -1) {
        // no distance between these two locations
        return NAN;
    }
    if (data->initial_tprops) {
        double tprop = data->initial_tprops[idx];
        return dist_from_tof(tprop - 0.5 * (data->delays[data->ds[a]] + data->delays[data->ds[b]]));
    } else {
        int peak_i = data->multilateration_peak_choices[idx];
        t_prop_peak_t *peak0;
        zarray_get_volatile(data->tprop_peaks[idx], peak_i, &peak0);
        return dist_from_tof(peak0->tprop - 0.5 * (data->delays[data->ds[a]] + data->delays[data->ds[b]]));
    }
}

void print_vec(int n, int per_row, double *vec)
{
    for (int i = 0; i < n; i++) {
        double v = vec[i];
        if (per_row == 1) {
            printf("\t");
        }
        if (fabs(v) < 1e-2 && v != 0) {
            printf("% 13.6e ", v);
        } else {
            printf("% 13.6f ", v);
        }
        if ((i % per_row) == (per_row - 1)) {
            printf("\n");
        }
    }
}

bool preconditioned_conjugate_gradient(int n, double *A, double *b, double epsilon, int kmax, double *x)
{
    // double zero_epsilon = 1e-100;

    double r[n];
    double z[n];
    double p[n];
    double aP[n];
    for (int i = 0; i < n; i++) {
        x[i] = 1.0;
        r[i] = 1.0;
    }

    bool done = false;
    int k = 0;

    // calculate residual
    for (int i = 0; i < n; i++) {
        double val = b[i];
        for (int j = 0; j < n; j++) {
            val -= A[i*n + j] * x[j];
        }
        r[i] = val;
    }

    double rDotR = 0.0;
    for (int i = 0; i < n; i++) {
        rDotR += r[i] * r[i];
    }

    double check = sqrt(rDotR);
    if (check < epsilon || k+1 >= kmax) {
        done = true;
    }

    double zDotR = 0.0;

    while (!done) {
        // this is the preconditioning step
        // we use the diagonal of A as 'Matrix' M
        // solve M*z = r for z.
        for (int i = 0; i < n; i++) {
            z[i] = r[i] / A[i*n+i];
            if (!isfinite(z[i])) {
                z[i] = r[i]; // just don't precondition this element then
            }
        }

        // on to the rest of the algorithm...
        double lastZDotR = zDotR;
        zDotR = 0.0;
        for (int i = 0; i < n; i++) {
            zDotR += z[i] * r[i];
        }

        if (k == 0) {
            memcpy(p, z, sizeof(z));
        } else {
            double delta = zDotR / lastZDotR;
            for (int i = 0; i < n; i++) {
                p[i] = z[i] + delta*p[i];
            }
        }
        double pAp = 0.0;
        for (int i = 0; i < n; i++) {
            double val = 0.0;
            for (int j = 0; j < n; j++) {
                val = val + A[i*n+j]*p[j];
            }
            aP[i] = val;
            pAp += p[i] * val;
        }

        double gamma = zDotR / pAp;

        for (int i = 0; i < n; i++) {
            x[i] += gamma * p[i];
        }

        for (int i = 0; i < n; i++) {
            r[i] -= gamma * aP[i];
        }

        // check for convergence again
        rDotR = 0.0;
        for (int i = 0; i < n; i++) {
            rDotR += r[i] * r[i];
        }

        check = sqrt(rDotR);
        if (check < epsilon || k+1 >= kmax) {
            done = true;
        }

        k++;
    }
    if (k >= kmax) {
        fprintf(stderr, "Failed to converge to solution of the linear system after %d iterations.\n", kmax);
        return false;
    } else {
        // printf("Converged to solution of linear system in %d iterations.\n", k);
    }
    // for (int i = 0; i < n; i++) {
    //     printf("%.3e ", x[i]);
    // }
    // printf("\n");

    return true;
}

double calc_neg_gradient_hessian_fval(int n, double (*f)(double*, void*), void* user, double *x, double *negGradient, double *hessian)
{
    // Since we need second derivatives, we have to use a larger value here.
    // the square root would be enough for the gradient, but so small that second derivative values would disappear
    // So we use the third root!
    double machineC = pow(1e-16, 1/3.0);

    double xPrime[n];

    double fVal = f(x, user);
    // printf("Current score: %.5f\n", fVal);

    double fNextVals[n];
    double fPrevVals[n];
    for (int i = 0; i < n; i++) {
        double theta = fabs(x[i])*machineC + machineC;
        memcpy(xPrime, x, sizeof(xPrime));
        xPrime[i] += theta;
        fNextVals[i] = f(xPrime, user);

        xPrime[i] -= theta * 2;
        fPrevVals[i] = f(xPrime, user);
    }

    // printf("Gradient: ");
    for (int i = 0; i < n; i++) {
        double theta = fabs(x[i])*machineC + machineC;
        // 2nd order centered finite difference for first derivative
        negGradient[i] = -(fNextVals[i] - fPrevVals[i]) / (2 * theta);
        // printf("%.3e ", -negGradient[i]);
    }
    // printf("\n");

    for (int i = 0; i < n; i++) {
        double theta1 = fabs(x[i])*machineC + machineC;
        // The Hessian is symmetric
        for (int j = i; j < n; j++) {
            if (i == j) {
                memcpy(xPrime, x, sizeof(xPrime));
                xPrime[i] -= theta1;
                double fPrime2 = f(xPrime, user);
                // 2nd order centered finite difference for second derivative
                hessian[i*n+j] = (fNextVals[i] + fPrime2 - 2*fVal) / (theta1 * theta1);
            } else {
                double theta2 = fabs(x[j])*machineC + machineC;
                memcpy(xPrime, x, sizeof(xPrime));
                xPrime[i] += theta1;
                xPrime[j] += theta2;
                double fBothNextVal = f(xPrime, user);
                // 2nd order forward finite difference for cross derivative
                hessian[i*n+j] = (fBothNextVal - fNextVals[i] - fNextVals[j] + fVal) / (theta1 * theta2);
                hessian[j*n+i] = hessian[i*n+j];
            }
        }
    }

    // printf("Hessian Matrix:\n");
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < n; j++) {
    //         printf("%8.1e ", hessian[i*n+j]);
    //     }
    //     printf("\n");
    // }

    return fVal;
}

double newton_optimize_base(int n, double (*f)(double*, void*), void* user, double *x, double epsilon, int kmax)
{
    double innerEpsilon = 1e-5;
    int inner_kmax = 250;
    double relaxation = 1.0;
    int relaxAttemptMax = 30;

    double newX[n];
    double negGradient[n];
    double hessian[n * n];
    double dx[n];
    int k = 0;
    double check = epsilon * 10;

    double fVal = 0.0;
    while (k < kmax && check > epsilon) {
        fVal = calc_neg_gradient_hessian_fval(n, f, user, x, negGradient, hessian);
        bool success = preconditioned_conjugate_gradient(n, hessian, negGradient, innerEpsilon, inner_kmax, dx);

        bool have_bad_value = false;
        for (int i = 0; i < n; i++) {
            if (!isfinite(dx[i])) {
                have_bad_value = true;
                break;
            }
        }
        if (!success || have_bad_value) {
            fprintf(stderr, "Warning: solution of Hx=-g had NaN values or did not converge. Terminating optimization early.\n");
            break;
        }

        double currentRelaxation = relaxation;
        int attemptNum = 0;
        double newFVal = 0.0;
        while (attemptNum < relaxAttemptMax) {
            for (int i = 0; i < n; i++) {
                newX[i] = x[i] + currentRelaxation*dx[i];
            }
            newFVal = f(newX, user);
            if (newFVal < fVal) {
                break;
            }
            currentRelaxation *= 0.5;
            attemptNum++;
        }
        if (newFVal >= fVal) {
            // printf("Relaxation down to %.3e did not improve value. Attempting individual element changes...\n", currentRelaxation);
            // return fVal;

            for (int i = 0; i < n; i++) {
                newX[i] = x[i];
            }
            double bestFVal = fVal;
            for (int i = 0; i < n; i++) {
                int attemptNum = 0;
                currentRelaxation = 1.0;
                while (attemptNum < relaxAttemptMax) {
                    newX[i] = x[i] + currentRelaxation*negGradient[i];
                    newFVal = f(newX, user);
                    if (newFVal < bestFVal) {
                        break;
                    }
                    currentRelaxation *= 0.5;
                    attemptNum++;
                }
                if (newFVal >= bestFVal) {
                    newX[i] = x[i]; // didn't improve, revert change
                } else {
                    bestFVal = newFVal;
                }
            }

            if (newFVal >= fVal) {
                // printf("Relaxation and individual element changes did not improve value. Convergence reached.\n");
                // printf("Converged to a minimum in %d iterations.\n", k + 1);
                return fVal;
            }
            newFVal = bestFVal;
            // printf("Using individual changes, fval: %.3e\n", newFVal);
        } else {
            // printf("Using relaxation: %.3e w/ fval: %.3e\n", currentRelaxation, newFVal);
        }

        memcpy(x, newX, sizeof(newX));

        check = fabs((fVal - newFVal) / fVal * 100);
        // printf("Score: %8.3f and percent change in score: %.3e\n", newFVal, check);
        k++;

        // printf("Iteration %d x values:\n", k);
        // for (int i = 0; i < n; i++) {
        //     printf("%.5f ", x[i]);
        // }
        // printf("\n");
    }

    if (k >= kmax) {
        fprintf(stderr, "Failed to converge to a minimum after %d iterations with check value %.4e\n", kmax, check);
    } else {
        // printf("Converged to a minimum in %d iterations.\n", k);
    }
    // for (int i = 0; i < n; i++) {
    //     printf("%.5f ", x[i]);
    // }
    // printf("\n");
    //
    // printf("Final score: %.5f\n", fVal);

    return fVal;
}

double newton_optimize(int n, double (*f)(double*, void*), void* user, double *x)
{
    double epsilon = 1e-3;
    int kmax = 70;
    return newton_optimize_base(n, f, user, x, epsilon, kmax);
}

void calc_errs(calib_data_t *data, double *coords, double *delays, double *errs)
{
    if (data->count_peak_choices) {
        printf("\n");
    }
    int below_3_stddev = 0;
    int above_3_stddev = 0;

    for (int i = 0; i < data->n_tprops; i++) {
        int a = data->as[i];
        int b = data->bs[i];

        double d = dist_between(data, coords, a, b);
        if (isinf(d)) {
            errs[i] = d; // propogate the infinity, this is a terrible solution!
            continue;
        }
        double dist_tprop = (d + range_bias_at(d)) / SPEED_OF_LIGHT +
                            0.5 * (delays[data->ds[a]] + delays[data->ds[b]]);

        if (data->densities) {
            measurement_density_t *density = &data->densities[i];
            // what is the likelihood of our dist_tprop in the density?
            double likelihood = 0;

            double max_x = density->min_x + (density->n - 1) * density->dx;
            if (dist_tprop < density->min_x || dist_tprop > max_x) {
                double diff = (dist_tprop - density->mean_x) / 1.3e-4;
                likelihood = -0.5 * diff * diff;
            } else {
                if ("quadratic") {
                    // use quadratic interpolation between points with precomputed coefficients
                    double loc = (dist_tprop - density->min_x) / density->dx;
                    int best_i = (int)(loc + 0.5);
                    double a = density->quadratic_as[best_i];
                    double b = density->quadratic_bs[best_i];
                    double c = density->quadratic_cs[best_i];
                    double p = a * dist_tprop * dist_tprop + b * dist_tprop + c;
                    likelihood = fastlog(p);

                    // double diff = (dist_tprop - data->initial_tprops[i]) / 1.3e-4;
                    // double l2 = -0.5 * diff * diff;
                    // // printf("%12.8e, %12.8e\n", likelihood, l2);
                    // likelihood = l2;

                    // if (i == 0) {
                    //     double diff = (dist_tprop - 0.5302234) / 1.3e-4;
                    //     double p2 = fastexp(-0.5 * diff * diff);
                    //     if ((p - p2) / p2 > 0.001) {
                    //         printf("\n%12.8e, %12.8e\n", p, p2);
                    //     }
                    // }
                } else {
                    // we do linear interpolation between points
                    double loc = (dist_tprop - density->min_x) / density->dx;
                    int low_i = (int)loc;
                    int high_i = low_i + 1;
                    if (high_i >= density->n) {
                        likelihood = fastlog(density->ys[low_i]);
                    } else {
                        double alpha = loc - low_i;
                        likelihood = fastlog(alpha * density->ys[high_i] + (1.0 - alpha) * density->ys[low_i]);
                    }
                }
            }

            // if (data->count_peak_choices && likelihood < -1) {
            //     printf("%d, %d: %.3f\n", a, b, likelihood);
            // }

            errs[i] = -likelihood;
        } else {
            zarray_t *tprop_peaks = data->tprop_peaks[i];
            int peaks_n = zarray_size(tprop_peaks);
            t_prop_peak_t *peaks = (t_prop_peak_t*)tprop_peaks->data;

            if (data->count_peak_choices) {
                // if (a == 6 && b == 1) {
                //     printf("\n12 and 81 dist_tprop: %.8e\n", dist_tprop);
                // }
                double diff = dist_tprop - peaks[0].tprop;
                double std_dev_off = diff / sqrt(peaks[0].variance);
                // printf("%d off by: %.3f\n", i, std_dev_off);
                if (std_dev_off > 3) {
                    above_3_stddev++;
                } else if (std_dev_off < -3) {
                    below_3_stddev++;
                }
            }

            int chosen_peak = 0;
            double max_likelihood = -DBL_MAX;
            // double product_likelihood = 0;
            // double sum_likelihood = 0;
            // double exp_sum_likelihood = 0;

            for (int j = 0; j < peaks_n; j++) {
                double diff = dist_tprop - peaks[j].tprop;
                double likelihood = -diff * diff / peaks[j].variance - 0.5 * peaks[j].log_variance + peaks[j].weight;
                // printf("diff %.6f variance: %.2e likelihood: %.4f\n", diff, peaks[j].variance, likelihood);
                if (likelihood > max_likelihood) {
                    if (data->count_peak_choices && j > 0) {
                        printf("\r");
                    }
                    max_likelihood = likelihood;
                    chosen_peak = j;
                }
                // product_likelihood += likelihood;
                // exp_sum_likelihood += exp(likelihood);
                // sum_likelihood = log(exp(likelihood) + exp(sum_likelihood));
            }
            // errs[i] = -product_likelihood;
            errs[i] = -max_likelihood;
            // errs[i] = -sum_likelihood;
            // errs[i] = log(exp_sum_likelihood);

            if (data->count_peak_choices) {
                if (chosen_peak == peaks_n - 1) {
                    data->chosen_peaks[15]++;
                } else {
                    data->chosen_peaks[chosen_peak]++;
                }
            }
        }
    }

    if (data->count_peak_choices && data->tprop_peaks) {
        printf("\n");
        for (int i = 0; i < 16; i++) {
            printf("%d ", data->chosen_peaks[i]);
        }

        printf("\n%.3f%% of time above by 3 standard deviations", 100.0 * above_3_stddev / data->n_tprops);
        printf("\n%.3f%% of time below by 3 standard deviations", 100.0 * below_3_stddev / data->n_tprops);
    }
}

// copies correct elements to fill up coords, but only places the correct pointer for delays
// since first three locations fix the coordinate frame,
// we fix the coordinates that are set to 0. (not even putting them in vars to be changed)
void extract_coords_delays(calib_data_t *data, double *vars, double *coords, double **delays)
{
    int n_locs = data->n_locs;
    coords[0] = 0;
    coords[1] = 0;
    coords[2] = vars[0];
    coords[3] = 0;
    memcpy(&coords[4], &vars[1], sizeof(double) * (2 * n_locs - 4));
    *delays = &vars[2 * n_locs - 3];
}

double norm_from_errors(calib_data_t *data, double *errs, double *delays, double *coords)
{
    double delays_normalizing = 0.0;
    double delay_sigma = data->delay_sigma;
    double inv_delay_sigma2 = 1.0 / (delay_sigma * delay_sigma);
    for (int i = 0; i < data->n_delays; i++) {
        double diff = delays[i] - data->expected_delay;
        delays_normalizing += 0.5 * diff * diff * inv_delay_sigma2;
    }

    double norm = 0;
    for (int i = 0; i < data->n_tprops; i++) {
        norm += errs[i];
    }

    // keeping the overall error values in the range of around 1 helps prevent numerical errors
    return (norm + delays_normalizing) / data->n_tprops;
}

double residual_at(double *vars, void *user)
{
    calib_data_t *data = user;
    double errs[data->n_tprops];
    double *delays;
    double *coords;

    int n_coords = data->n_locs * 2;
    double coords_data[n_coords];

    if (data->opt_delays && data->opt_coords) {
        extract_coords_delays(data, vars, coords_data, &delays);
        coords = coords_data;
        calc_errs(data, coords, delays, errs);
    } else if (data->opt_coords) {
        extract_coords_delays(data, vars, coords_data, &delays);
        coords = coords_data;
        delays = data->delays;
        calc_errs(data, coords, data->delays, errs);
    } else if (data->opt_delays) {
        delays = vars;
        coords = data->coords;
        calc_errs(data, data->coords, vars, errs);
    } else {
        return 0;
    }

    return norm_from_errors(data, errs, delays, coords);
}

// false on failure due to a singular matrix, points that poorly define the axis
// will likely only fail then with 3 or 4 points chosen in the mask.
// https://people.eecs.berkeley.edu/~prabal/projects/cs294-1/multilateration.pdf
bool multilaterate(int n_locs, double *coords, bool *mask, double *rs, double *x, double *y)
{
    // exactly one set of coordinates are used as the offset for the entire multilateration
    // so it is important to have a good low-error choice. We just use the closest value
    int n_vals = 0;

    int min_rs2_idx = -1;
    double min_rs2 = DBL_MAX;
    for (int i = 0; i < n_locs; i++) {
        if (!mask[i]) {
            continue;
        }
        n_vals++;
        rs[i] = rs[i] * rs[i];
        if (rs[i] < min_rs2) {
            min_rs2 = rs[i];
            min_rs2_idx = i;
        }
    }
    int comp_idx = min_rs2_idx;

    double min_x = coords[comp_idx*2];
    double min_y = coords[comp_idx*2+1];

    double b_const = min_x * min_x + min_y * min_y - rs[comp_idx];
    mask[comp_idx] = false;
    n_vals--;

    double A[n_vals * 2];
    double b[n_vals];

    int val_i = 0;
    for (int i = 0; i < n_locs; i++) {
        if (!mask[i]) {
            continue;
        }

        double x_i = coords[i*2];
        double y_i = coords[i*2+1];
        A[val_i*2] = min_x - x_i;
        A[val_i*2+1] = min_y - y_i;

        b[val_i] = 0.5 * (rs[i] - x_i * x_i - y_i * y_i + b_const);
        val_i++;
    }

    double ATA[4]; // 2x2 A^T*A
    double ATb[2]; // 2x1 A^T*b
    for (int a = 0; a < 2; a++) {
        for (int b = 0; b <= a; b++) {
            double ata_sum = 0;
            for (int i = 0; i < n_vals; i++) {
                ata_sum += A[i * 2 + a] * A[i * 2 + b];
            }
            ATA[a * 2 + b] = ata_sum;
            ATA[b * 2 + a] = ata_sum;
        }
        double atb_sum = 0;
        for (int i = 0; i < n_vals; i++) {
            atb_sum += A[i * 2 + a] * b[i];
        }
        ATb[a] = atb_sum;
    }

    // Gaussian elimination of ATA * x = Atb
    double xy[2];
    double factor = -ATA[2] / ATA[0];
    ATA[3] += factor * ATA[1];
    ATb[1] += factor * ATb[0];
    xy[1] = ATb[1] / ATA[3];
    xy[0] = (ATb[0] - ATA[1] * xy[1]) / ATA[0];

    if (!isfinite(xy[0]) || !isfinite(xy[1])) {
        mask[comp_idx] = true; // set back to how it was given to us
        return false;
    }

    *x = xy[0];
    *y = xy[1];

    mask[comp_idx] = true; // set back to how it was given to us
    return true;
}

void multilaterate_with(calib_data_t *data, int i, double *x, double *y)
{
    int n_locs = data->n_locs;
    double rs[n_locs];
    bool mask[n_locs];
    for (int j = 0; j < data->n_locs; j++) {
        // only valid to use actual different beacons
        // if two indices use the same delay time, they are the same beacon
        mask[j] = false;
        if (data->ds[j] != data->ds[i]) {
            double d = dist_for(data, j, i);
            if (!isnan(d)) {
                rs[j] = d;
                mask[j] = true;
            }
        }
    }

    multilaterate(n_locs, data->coords, mask, rs, x, y);
}

void iterate_multilateration(calib_data_t *data)
{
    for (int i = 3; i < data->n_locs; i++) {
        multilaterate_with(data, i, &data->coords[i*2], &data->coords[i*2+1]);
    }
}

// 0 if not robust, otherwise the test value used
double is_robust_triangle(double *coords, int *is)
{
    double min_length = DBL_MAX;
    double min_angle = DBL_MAX;
    double atans[3];
    for (int i = 0; i < 3; i++) {
        int this_i = is[i];
        int next_i = is[(i + 1) % 3];
        double dx = coords[this_i * 2] - coords[next_i * 2];
        double dy = coords[this_i * 2 + 1] - coords[next_i * 2 + 1];
        atans[i] = atan2(dy, dx);

        double l = sqrt(dx * dx + dy * dy);
        if (l < min_length) {
            min_length = l;
        }
    }

    for (int i = 0; i < 3; i++) {
        int next_i = (i + 1) % 3;
        double a = fabs(atans[i] - atans[next_i]);
        if (a > M_PI) {
            a = 2 * M_PI - a;
        }
        if (a < min_angle) {
            min_angle = a;
        }
    }

    double sin_min_angle = sin(min_angle);
    double test_val = min_length * sin_min_angle * sin_min_angle;
    // printf("test_val: %.3f\n", test_val);
    if (test_val < 0.01) { //if (test_val < 0.001) {
        return 0;
    }
    return test_val;
}

// 0 if not robust, otherwise the minimum of test values used
double is_robust_quad(double *coords)
{
    double min_val = DBL_MAX;
    for (int i = 0; i < 4; i++) {
        int is[3] = { 0 };
        for (int j = 0; j < 3; j++) {
            is[j] = (i + j) % 4;
        }
        double robust = is_robust_triangle(coords, is);
        if (robust == 0) {
            return 0;
        }
        if (robust < min_val) {
            min_val = robust;
        }
    }
    return min_val;
}

// false if it did not exceed because a distance did not exist between locations in `is`
bool quad_triangulation(calib_data_t *data, double *coords, int *is)
{
    // by definition: second position marks x direction
    double x2 = dist_for(data, is[0], is[1]);

    // by definition: first position
    coords[0] = 0;
    coords[1] = 0;

    coords[2] = x2;
    coords[3] = 0;

    // by definition: third position is on positive y side of line between 1st and 2nd locations
    double a = dist_for(data, is[0], is[2]);
    double b = dist_for(data, is[1], is[2]);

    if (isnan(a) || isnan(b) || isnan(x2)) {
        // fprintf(stderr, "Fatal error in quad_triangulation: distances did not exist between base locations\n");
        // exit(1);
        return false;
    }

    double x3 = (a*a - b*b + x2*x2) / (2 * x2);
    double y3 = sqrt(fabs(a*a - x3*x3));

    coords[4] = x3;
    coords[5] = y3;

    double d1 = dist_for(data, is[0], is[3]);
    double d2 = dist_for(data, is[1], is[3]);
    double d3 = dist_for(data, is[2], is[3]);

    if (isnan(d1) || isnan(d2) || isnan(d3)) {
        return false;
    }

    double x4 = (x2*x2 + d1*d1 - d2*d2) / (2 * x2);
    double y4 = (x3*x3 + y3*y3 - 2*x3*x4 + d1*d1 - d3*d3) / (2 * y3);
    coords[6] = x4;
    coords[7] = y4;

    return true;
}

double find_robust_quad_for(calib_data_t *data, int i, int n_solved_locs, int *is)
{
    int n = data->n_locs;
    bool mask[n];
    memset(mask, 0, sizeof(mask));
    double rs[n];

    double max_score = 0;
    double max_score_x;
    double max_score_y;

    for (int i1 = 0; i1 < n_solved_locs; i1++) {
        rs[is[i1]] = dist_for(data, is[i1], i);
        if (isnan(rs[is[i1]])) {
            continue;
        }
        mask[is[i1]] = true;
        for (int i2 = i1 + 1; i2 < n_solved_locs; i2++) {
            rs[is[i2]] = dist_for(data, is[i2], i);
            if (isnan(rs[is[i2]])) {
                continue;
            }
            mask[is[i2]] = true;
            for (int i3 = i2 + 1; i3 < n_solved_locs; i3++) {
                rs[is[i3]] = dist_for(data, is[i3], i);
                if (isnan(rs[is[i3]])) {
                    continue;
                }
                mask[is[i3]] = true;

                // At this point we have found 3 already solved points
                // And we can try to triangulate our new one and check if it is robust
                double new_x;
                double new_y;
                bool success = multilaterate(n, data->coords, mask, rs, &new_x, &new_y);
                if (success) {
                    double quad_coords[8];
                    quad_coords[0] = data->coords[is[i1] * 2];
                    quad_coords[1] = data->coords[is[i1] * 2 + 1];
                    quad_coords[2] = data->coords[is[i2] * 2];
                    quad_coords[3] = data->coords[is[i2] * 2 + 1];
                    quad_coords[4] = data->coords[is[i3] * 2];
                    quad_coords[5] = data->coords[is[i3] * 2 + 1];
                    quad_coords[6] = new_x;
                    quad_coords[7] = new_y;
                    double score = is_robust_quad(quad_coords);
                    // printf("%d: score: %.4f\n", i, score);
                    if (score > max_score) {
                        max_score = score;
                        max_score_x = new_x;
                        max_score_y = new_y;
                        goto outside;
                        // printf("%d: score: %.4f\n", i, max_score);
                    }
                }

                mask[is[i3]] = false;
            }
            mask[is[i2]] = false;
        }
        mask[is[i1]] = false;
    }
    outside:
    if (max_score > 0) {
        data->coords[i * 2] = max_score_x;
        data->coords[i * 2 + 1] = max_score_y;
    }
    return max_score;
}

bool find_first_robust_quad(calib_data_t *data, int *is, double *quad_coords)
{
    int n = data->n_locs;
    for (int i = 0; i < n; i++) {
        is[0] = i;
        for (int j = i + 1; j < n; j++) {
            is[1] = j;
            for (int k = j + 1; k < n; k++) {
                is[2] = k;
                for (int l = k + 1; l < n; l++) {
                    is[3] = l;
                    bool success = quad_triangulation(data, quad_coords, is);
                    if (success && is_robust_quad(quad_coords)) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

double inner_reset_multilateration(calib_data_t *data, bool *all_success)
{
    *all_success = true;

    int n = data->n_locs;

    double quad_coords[8];
    bool solved_locs[n];
    memset(solved_locs, 0, sizeof(solved_locs));
    int n_solved_locs = 0;

    // find the first robust quad we can
    int is[n];
    bool found_robust = find_first_robust_quad(data, is, quad_coords);
    if (!found_robust) {
        fprintf(stderr, "Error in reset_multilateration: no robust quads found\n");
        exit(1);
    }

    for (int i = 0; i < 4; i++) {
        data->coords[is[i] * 2] = quad_coords[i * 2];
        data->coords[is[i] * 2 + 1] = quad_coords[i * 2 + 1];
        solved_locs[is[i]] = true;
        n_solved_locs++;
    }

    // now we want to incrementally add single locations, but only if they make a robust quad
    // with some three other already solved points
    for (int i = 0; i < n; i++) {
        if (solved_locs[i]) {
            continue;
        }

        bool success = find_robust_quad_for(data, i, n_solved_locs, is);
        if (!success) {
            *all_success = false;
            // fprintf(stderr, "Warning in reset_multilateration: failed to find a robust quad; solution may be suboptimal\n");
        }

        solved_locs[i] = true;
        is[n_solved_locs] = i;
        n_solved_locs++;
    }

    double errs[data->n_tprops];
    calc_errs(data, data->coords, data->delays, errs);
    return norm_from_errors(data, errs, data->delays, data->coords);
}

// finding the best initial guess mostly helps speed up the non-linear solver
// whether it affects the final result mainly depends on the basic of convergence
// of the objective function (so it shouldn't with a well-formulated function)
double reset_multilateration(calib_data_t *data)
{
    int n_tprops = data->n_tprops;
    int best_peaks[n_tprops];
    memset(best_peaks, 0, sizeof(best_peaks));

    int peak_choices[n_tprops];
    memset(peak_choices, 0, sizeof(peak_choices));
    data->multilateration_peak_choices = peak_choices;

    bool success = false;
    double residual = inner_reset_multilateration(data, &success);
    if (data->densities) {
        return residual;
    }
    // return residual;
    // printf("initial residual %.6e\n", residual);

    for (int round_i = 0; round_i < 2; round_i++) {
        for (int i = 0; i < n_tprops; i++) {
            zarray_t *peaks = data->tprop_peaks[i];
            // t_prop_peak_t *peak_data = (t_prop_peak_t*)peaks->data;
            int peak_n = zarray_size(peaks);
            if (peak_n <= 1) {
                continue;
            }

            for (int j = 0; j < peak_n; j++) {
                // if (fabs(peak_data[j].tprop - peak_data[best_peaks[i]].tprop) <= 1e-5) {
                //     continue;
                // }
                peak_choices[i] = j;
                // printf("%e\n",    fabs(peak_data[j] - peak_data[0]));
                bool new_success = false;
                double new_residual = inner_reset_multilateration(data, &new_success);
                // printf("i %2d j %d new_residual %.6e\n", i, j, new_residual);

                if ((!success && new_success) || (new_success && (new_residual < residual))) {
                    // return new_residual;
                    residual = new_residual;
                    success = new_success;
                    best_peaks[i] = j;
                }
            }
            peak_choices[i] = best_peaks[i];
        }
    }

    if (!success) {
        fprintf(stderr, "Warning in reset_multilateration: failed to find a robust quad; solution may be suboptimal\n");
    }

    // printf("final_residual %.6e\n", residual);
    return residual;
}

double perform_optimize(calib_data_t *data)
{
    double residual;

    int n_locs = data->n_locs;
    int n_delays = data->n_delays;
    // big enough for all cases
    int n_coords = n_locs * 2;
    double passing_vars[n_coords - 3 + n_delays];

    double *coords = data->coords;
    double *delays = data->delays;

    if (data->opt_delays && !data->opt_coords) {
        residual = newton_optimize_base(n_delays, residual_at, data, delays, 1e-5, 1000);
        return residual;
    } else if (!data->opt_delays && !data->opt_coords) {
        // do nothing
        return NAN;
    }

    passing_vars[0] = coords[2];
    memcpy(&passing_vars[1], &coords[4], sizeof(double) * (2 * n_locs - 4));
    if (data->opt_delays) {
        memcpy(&passing_vars[2 * n_locs - 3], delays, sizeof(double) * n_delays);
        residual = newton_optimize_base(2 * n_locs - 3 + n_delays, residual_at, data, passing_vars, 1e-5, 1000);
        memcpy(delays, &passing_vars[2 * n_locs - 3], sizeof(double) * n_delays);
    } else {
        residual = newton_optimize_base(2 * n_locs - 3, residual_at, data, passing_vars, 1e-5, 1000);
    }
    coords[2] = passing_vars[0];
    memcpy(&coords[4], &passing_vars[1], sizeof(double) * (2 * n_locs - 4));

    return residual;
}

calib_data_t *beacon_autocal_create(int n_locs, int n_delays, int n_tprops,
                                    int z_plane_delay_idx, double z_plane_2d_offset,
                                    measurement_density_t *densities, int *as, int *bs, int *ds,
                                    double delay_mean, double delay_sigma, bool recalibrate_delays)
{
    bool opt_delays = n_locs >= 7 && recalibrate_delays;
    calib_data_t data = { .n_locs = n_locs, .n_delays = n_delays, .n_tprops = n_tprops,
                          .densities = densities, .as = as, .bs = bs, .ds = ds,
                          .opt_coords = true, .opt_delays = opt_delays,
                          .expected_delay = delay_mean, .delay_sigma = delay_sigma,
                          .z_plane_delay_idx = z_plane_delay_idx, .z_plane_2d_offset = z_plane_2d_offset };

    // precompute quadratic fit for range bias calculation
    initialize_range_bias_quadratic_fit();

    calib_data_t *ret_data = malloc(sizeof(data));
    memcpy(ret_data, &data, sizeof(data));

    return ret_data;
}

void beacon_autocal_destroy(calib_data_t *data)
{
    free(data);
}

double beacon_autocal_trilaterate(calib_data_t *data, zarray_t **tprop_peaks, double *initial_tprops, double *delays, double *coords)
{
    data->tprop_peaks = tprop_peaks;
    data->initial_tprops = initial_tprops;
    data->delays = delays;
    data->coords = coords;

    if (data->opt_delays) {
        for (int i = 0; i < data->n_delays; i++) {
            delays[i] = data->expected_delay;
        }
    }

    double residual;
    if (data->opt_coords) {
        residual = reset_multilateration(data);
    } else {
        double errs[data->n_tprops];
        calc_errs(data, data->coords, data->delays, errs);
        residual = norm_from_errors(data, errs, delays, coords);
    }

    return residual;
}

double beacon_autocalibrate(calib_data_t *data, zarray_t **tprop_peaks, double *initial_tprops, double *delays, double *coords, bool verbose)
{
    data->tprop_peaks = tprop_peaks;
    data->initial_tprops = initial_tprops;
    data->delays = delays;
    data->coords = coords;

    if (data->opt_delays) {
        for (int i = 0; i < data->n_delays; i++) {
            delays[i] = data->expected_delay;
        }
    }

    // precompute log variances
    if (tprop_peaks) {
        for (int i = 0; i < data->n_tprops; i++) {
            zarray_t *peaks = tprop_peaks[i];
            int peaks_n = zarray_size(peaks);
            t_prop_peak_t *peak_data = (t_prop_peak_t*)peaks->data;
            for (int j = 0; j < peaks_n; j++) {
                peak_data[j].log_variance = log(peak_data[j].variance);
            }
        }
    }

    double residual;
    if (data->opt_coords) {
        residual = reset_multilateration(data);
    }  else {
        double errs[data->n_tprops];
        calc_errs(data, data->coords, data->delays, errs);
        residual = norm_from_errors(data, errs, delays, coords);
    }
    if (verbose) {
        printf("Starting residual: %.4f\n", residual);
    }

    // Full non-linear optimization
    residual = perform_optimize(data);

    printf("\r%11s %11.6f\r", "FINAL", residual);
    fflush(stdout);

    if (verbose) {
        printf("\n");
        // print_vec(n_delays, n_delays, delays);
        // printf("Errs:\n");
        // print_vec(n_tprops, errs);
        printf("Coords:\n");
        print_vec(data->n_locs * 2, 2, coords);
    }

    if (verbose) {
        data->count_peak_choices = true;
        double errs[data->n_tprops];
        calc_errs(data, data->coords, data->delays, errs);
        data->count_peak_choices = false;
    }

    return residual;
}
