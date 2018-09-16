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

#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>
#include <sys/select.h>

#include "beacon_calibration.h"
#include "common/getopt.h"
#include "common/time_util.h"

#define PEAK_DELTA 5e-6
#define MERGE_EPSILON 5e-6 // ~0.15cm

#define MAX_BEACONS 255

#define SIGMA 1.3e-4
#define INV_SIGMA (1.0 / SIGMA)
#define NARROW_SIGMA 1.6e-5
#define INV_NARROW_SIGMA (1.0 / NARROW_SIGMA)

#define WEIGHT_BROAD -40
#define VARIANCE_BROAD 1e3

#define DEFAULT_ANTENNA_DELAY 0.516

typedef struct {
    getopt_t *gopt;

    // exit with current measurements regardless of when we otherwise would. helpful when using a dataset file.
    bool take_no_more_measurements;
    // debugging with a log file that doesn't have all the measurements
    bool no_all_measurements;

    // Uses the variance of an entire set of pairwise measurements instead of the variance of each individual
    // peak as calulated by finding the full width at half max.
    bool use_whole_measurement_variance;

    // instead of finding peaks in the data, we estimate the probability density of the measurements
    // and use that directly instead, using the peaks in it only as initial conditions
    bool use_measurement_density;

    // number of full rounds (each containing all those perturbations for choosing initialization)
    int n_full_rounds;

    // number of randomly sampled initializations for measurement density based inits.
    int n_perturbation_rounds;

    // enables loading of additional locations where the robot has moved after the initial one
    bool use_dance_locations;

    // prior on the antenna delays
    double delay_mean;
    double delay_sigma;

    bool evaluate_by_measurements;
    bool evaluate_by_full_rounds;
    bool evaluate_by_pert_rounds;
    bool evaluate_by_dancing;
    bool evaluate_by_delay_mean;
    bool evaluate_by_delay_sigma;

    // evaluation parameters
    int eval_dataset;
    bool first_peak;
    bool max_peak;
    bool triangulation_only;
    bool fixed_delays;
    int eval_run_i;
    double ground_truth_xs[32];
    double ground_truth_ys[32];
    double eval_min_x;
    double eval_max_x;
    double eval_min_y;
    double eval_max_y;

    // maps beacon id number to location index in xs and ys
    // -1 if the beacon is not present
    int beacon_loc_i[MAX_BEACONS];

    bool try_alternate_tprops;
    double robot_z_offset; // in 2d case, optional z-offset between robot and all other beacons

    bool silent_output;

    int n_locs;
    int n_locs_alloc;
    int n_beacons; // equivalent to n_delays
    int n_tprops;
    int n_tprops_alloc;

    int *as; // maps tprop index to "a" location index in xs and ys, n_tprops used, n_tprops_alloc total
    int *bs; // maps tprop index to "b" location index in xs and ys, n_tprops used, n_tprops_alloc total
    int *ds; // maps location index to antenna delay index, n_locs used, n_locs_alloc total
    int *ids; // maps location index to beacon id, n_locs used, n_locs_alloc total

    double *delays; // n_beacons used, n_locs_alloc total
    double *xs; // n_locs used, n_locs_alloc total
    double *ys; // n_locs used, n_locs_alloc total

    bool recalibrate_delays;

    int self_id;
    all_measurements_t *all_measurements; // n_tprops used, n_tprops_alloc total
    // probability density calculated from all_measurements
    measurement_density_t *densities; // n_tprops used, n_tprops_alloc total

    int measurements_removed;
    int measure_epsilon_times_reached;
    zarray_t **all_t_prop_peaks; // n_tprops used, n_tprops_alloc total, zarray of t_prop_peak_t

    bool estimation_ready;
    bool estimation_complete;
} beacon_localize_state_t;

// Adds an entry for a new propagation time in the calibration system between locations a and b
// now we have another measurement to make and use for solving the calibration
void add_new_tprop(beacon_localize_state_t *state, int a, int b)
{
    if (state->n_tprops >= state->n_tprops_alloc) {
        int new_alloc = state->n_tprops_alloc * 2;
        new_alloc = new_alloc < 32 ? 32 : new_alloc;

        state->as = realloc(state->as, sizeof(int) * new_alloc);
        state->bs = realloc(state->bs, sizeof(int) * new_alloc);
        state->all_t_prop_peaks = realloc(state->all_t_prop_peaks, sizeof(zarray_t*) * new_alloc);
        state->all_measurements = realloc(state->all_measurements, sizeof(all_measurements_t) * new_alloc);
        state->densities = realloc(state->densities, sizeof(measurement_density_t) * new_alloc);

        state->n_tprops_alloc = new_alloc;
    }

    int idx = state->n_tprops;
    state->n_tprops++;

    state->as[idx] = a;
    state->bs[idx] = b;
    state->all_t_prop_peaks[idx] = zarray_create(sizeof(t_prop_peak_t));
    state->all_measurements[idx].n = 0;
    state->all_measurements[idx].a = a;
    state->all_measurements[idx].b = b;
    for (int i = 0; i < MAX_MEASUREMENTS; i++) {
        state->all_measurements[idx].tprops[i] = 1337.42;
    }
    state->densities[idx].n = 0;
    state->densities[idx].ys = NULL;
    state->densities[idx].quadratic_as = NULL;
    state->densities[idx].quadratic_bs = NULL;
    state->densities[idx].quadratic_cs = NULL;
}

// Adds an entry for a new location in the calibration system
// now we have another location to measure and solve for later
void add_new_loc(beacon_localize_state_t *state, int id)
{
    if (state->n_locs >= state->n_locs_alloc) {
        int new_alloc = state->n_locs_alloc * 2;
        new_alloc = new_alloc < 32 ? 32 : new_alloc;

        state->ds = realloc(state->ds, sizeof(int) * new_alloc);
        state->ids = realloc(state->ids, sizeof(int) * new_alloc);
        state->delays = realloc(state->delays, sizeof(double) * new_alloc);
        state->xs = realloc(state->xs, sizeof(double) * new_alloc);
        state->ys = realloc(state->ys, sizeof(double) * new_alloc);

        state->n_locs_alloc = new_alloc;
    }

    int idx = state->n_locs;
    state->n_locs++;

    // have we used this specific id before? If not, add a new beacon/delay time
    if (state->beacon_loc_i[id] == -1) {
        int delay_idx = state->n_beacons;
        state->n_beacons++;

        state->delays[delay_idx] = 0;
        state->ds[idx] = delay_idx;

        // prep to make measurements to all previously existing beacons
        for (int b = 0; b < idx; b++) {
            add_new_tprop(state, idx, b);
        }
    } else {
        // get delay index from the previous instance of this beacon id
        state->ds[idx] = state->ds[state->beacon_loc_i[id]];

        // prep to make measurements to all other beacons
        for (int i = 0; i < MAX_BEACONS; i++) {
            int b = state->beacon_loc_i[i];
            if (i != id && b != -1) {
                add_new_tprop(state, idx, b);
            }
        }
    }

    state->beacon_loc_i[id] = idx;
    state->ids[idx] = id;

    state->xs[idx] = 0;
    state->ys[idx] = 0;
}

int double_cmp(const void *a, const void *b)
{
    if (*(double*)a <  *(double*)b) {
        return -1;
    } else if (*(double*)a == *(double*)b){
        return 0;
    } else {
        return 1;
    }
}

double sort_and_median(all_measurements_t *all)
{
    qsort(all->tprops, all->n, sizeof(double), double_cmp);
    return all->tprops[all->n / 2];
}

double eval_gauss_peak(int n, double *xs, double mu, double sigma, double inv_sigma)
{
    double low_bound = -4 * sigma + mu;
    double high_bound = 4 * sigma + mu;
    double val = 0;
    for (int i = 0; i < n; i++) {
        if (xs[i] < low_bound || xs[i] > high_bound) {
            continue;
        }

        double a = (xs[i] - mu) * inv_sigma;
        val += fastexp(-0.5 * a * a);
    }
    return val / n;
}

double eval_gauss_peak_all(int n, double x, double *mus, double inv_sigma)
{
    double val = 0;
    for (int i = 0; i < n; i++) {
        double a = (x - mus[i]) * inv_sigma;
        val += fastexp(-0.5 * a * a);
    }
    return val;
}

typedef struct
{
    int n;
    double *xs;
    double *ys;
    double max_val;
    double mu;
} sigma_fit_t;

double sigma_fit_f(double *sigma, void *user)
{
    sigma_fit_t *data = (sigma_fit_t*)user;
    double sig = *sigma;
    double inv_sig = 1.0 / sig;

    double err_sum = 0;
    // double c = (1.0 / sqrt(2 * M_PI)) * inv_sig;
    // double normalizing_c = max_val / c;

    for (int i = 0; i < data->n; i++) {
        double x = data->xs[i];
        double y = data->ys[i];

        double a = (x - data->mu) * inv_sig;
        double expected_val = data->max_val * fastexp(-0.5 * a * a);

        double diff = expected_val - y;
        err_sum += diff * diff;
    }

    return err_sum;
}

double find_weight_for_gaussian(int n, double *xs, double mu, double variance)
{
    double std = sqrt(variance);
    double low = mu - std;
    double high = mu + std;

    int count = 0;
    for (int i = 0; i < n; i++) {
        if (xs[i] >= low && xs[i] <= high) {
            count++;
        }
    }

    // printf("mu: %.8e variance: %.8e count: %d\n", mu, variance, count);
    return log((double)count / n);
}

double fit_gaussian_variance(int n, double *xs, double mu)
{
    int smoothed_n = 400; // higher values cause problems... why?
    double smoothed_xs[smoothed_n];
    double smoothed_ys[smoothed_n];

    double max_val = eval_gauss_peak(n, xs, mu, NARROW_SIGMA, INV_NARROW_SIGMA);
    double min_x = DBL_MAX;
    double max_x = -DBL_MAX;

    for (int i = 0; i < n; i++) {
        if (xs[i] < min_x) {
            min_x = xs[i];
        }
        if (xs[i] > max_x) {
            max_x = xs[i];
        }
    }

    for (int i = 0; i < smoothed_n; i++) {
        double x = min_x + (max_x - min_x) / smoothed_n * i;
        smoothed_xs[i] = x;
        smoothed_ys[i] = eval_gauss_peak(n, xs, x, NARROW_SIGMA, INV_NARROW_SIGMA);
    }

    sigma_fit_t fit_data = { .n = smoothed_n, .xs = smoothed_xs, .ys = smoothed_ys, .max_val = max_val, .mu = mu };
    double new_sigma = SIGMA;
    newton_optimize(1, sigma_fit_f, &fit_data, &new_sigma);

    return new_sigma * new_sigma;
}

void inner_find_peaks(int n, double *xs, double start_mu, double sigma, zarray_t *peaks, t_prop_peak_t *peak_max)
{
    // Peaks are returned either in the given peaks zarray, OR, just the median-based in peak_median
    // depending on whether peaks in null or not

    double delta = PEAK_DELTA;
    double big_delta = delta * 32;

    double inv_sigma = 1.0 / sigma;
    double thresh = 0.01;

    // peak max just finds the max around start_mu
    // otherwise, we scan the whole set of xs
    int steps_n = 1;
    if (!peak_max) {
        double min_x = DBL_MAX;
        double max_x = -DBL_MAX;
        for (int i = 0; i < n; i++) {
            double v = xs[i];
            if (v < min_x) {
                min_x = v;
            }
            if (v > max_x) {
                max_x = v;
            }
        }

        steps_n = (int)((max_x - min_x) / big_delta + 0.5) + 1;
        start_mu = max_x;
        // printf("steps: %d, n: %d\n", steps_n, n);
    }

    for (int i = 0; i < steps_n; i++) {
        double mu = start_mu - big_delta * i;
        double best_mu = mu;
        double best_fval = eval_gauss_peak(n, xs, mu, sigma, inv_sigma);
        // printf("%d: mu %.13lf fval %.13lf\n", i, mu, best_fval);
        bool reached_peak_low = false;
        bool reached_peak_high = false;

        int inner_n = (int)(big_delta / delta + 0.5);
        for (int j = 0; j < inner_n; j++) {
            mu -= delta;
            // printf("\t%d: mu %.13lf\n", j, mu);
            double new_fval = eval_gauss_peak(n, xs, mu, sigma, inv_sigma);
            if (new_fval > best_fval) {
                best_fval = new_fval;
                best_mu = mu;
            } else {
                reached_peak_low = true;
                // printf("low at %d mu: %f\n", j, mu);
                break;
            }
        }

        mu = start_mu - big_delta * i;
        for (int j = 0; j < inner_n; j++) {
            mu += delta;
            double new_fval = eval_gauss_peak(n, xs, mu, sigma, inv_sigma);
            if (new_fval > best_fval) {
                best_fval = new_fval;
                best_mu = mu;
            } else {
                reached_peak_high = true;
                break;
            }
        }

        if (peak_max) {
            peak_max->tprop = best_mu;
            peak_max->variance = fit_gaussian_variance(n, xs, best_mu);
            peak_max->weight = find_weight_for_gaussian(n, xs, best_mu, peak_max->variance);
            return;
        }

        if (reached_peak_low && reached_peak_high && best_fval > thresh) {
            double variance = fit_gaussian_variance(n, xs, best_mu);
            double weight = find_weight_for_gaussian(n, xs, best_mu, variance);
            t_prop_peak_t peak = { .tprop = best_mu, .variance = variance, .weight = weight };
            zarray_add(peaks, &peak);
        }
    }

    if (zarray_size(peaks) == 0) {
        exit(1);
    }
}

void find_peaks(int n, double *xs, double variance, zarray_t *peaks, int *modes_n)
{
    // Really dumb search and gradient ascent to find best gaussian mean mu
    // look for local minima to consider as peaks
    // Although most peaks are based on the empirically found expected standard deviation,
    // a peak is also included that is the essentially the maximum location
    // which may be useful when a unimodal distribution is skewed.
    // The number of modes is also returned, which may be less than or equal to the number of peaks.
    // This is a lot of heuristic guessing, so it can definitely be improved at some point...
    zarray_clear(peaks);

    inner_find_peaks(n, xs, 0, SIGMA, peaks, NULL);

    t_prop_peak_t* data = (t_prop_peak_t*)peaks->data;
    zarray_sort(peaks, double_cmp);

    data = (t_prop_peak_t*)peaks->data;
    int modes = zarray_size(peaks);

    // remove effectively identical peaks
    for (int i = 0; i < zarray_size(peaks) - 1; i++) {
        double diff = fabs(data[i].tprop - data[i + 1].tprop);
        if (diff <= SIGMA / 2) {
            modes--;
        }
        if (diff <= MERGE_EPSILON) {
            zarray_remove_index(peaks, i + 1, false);
            i--;
        }
    }

    int base_peaks_n = zarray_size(peaks);
    if (base_peaks_n > 0) {
        double avg_tprop = 0;
        // for (int i = 0; i < n; i++) {
        //     avg_tprop += xs[i];
        // }
        // avg_tprop /= n;
        for (int i = 0; i < zarray_size(peaks); i++) {
            avg_tprop += data[i].tprop;
        }
        avg_tprop /= zarray_size(peaks);
        t_prop_peak_t avg_peak = { .tprop = avg_tprop, .variance = variance * VARIANCE_BROAD, .weight = WEIGHT_BROAD }; // 2e-3, or -4 or -5
        zarray_add(peaks, &avg_peak);
    }

    data = (t_prop_peak_t*)peaks->data;

    zarray_sort(peaks, double_cmp);

    if (modes_n) {
        *modes_n = modes;
    }
}

void calculate_qudratic_density_fit(measurement_density_t *density)
{
    int n = density->n;

    density->quadratic_as = realloc(density->quadratic_as, sizeof(double) * n);
    density->quadratic_bs = realloc(density->quadratic_bs, sizeof(double) * n);
    density->quadratic_cs = realloc(density->quadratic_cs, sizeof(double) * n);

    // double inv_4_dx2 = 1.0 / (4.0 * density->dx * density->dx);
    double inv_2_dx = 1.0 / (2.0 * density->dx);
    // double inv_dx = 1.0 / (density->dx);
    double inv_dx2 = 1.0 / (density->dx * density->dx);

    for (int i = 1; i < n - 1; i++) {
        double x2 = density->min_x + density->dx * i;
        double y1 = density->ys[i - 1];
        double y2 = density->ys[i];
        double y3 = density->ys[i + 1];
        // double a = (y3 - 2 * y2 + y1) * inv_4_dx2;
        double a = (0.5 * (y3 + y1) - y2) * inv_dx2;
        double b = (y3 - y1) * inv_2_dx - 2.0 * a * x2;
        double c = y2 - a * x2 * x2 - b * x2;

        density->quadratic_as[i] = a;
        density->quadratic_bs[i] = b;
        density->quadratic_cs[i] = c;

        // double y2_2 = a * x2 * x2 + b * x2 + c;
        // printf("%.8f\n", y2_2);
    }

    density->quadratic_as[0] = density->quadratic_as[1];
    density->quadratic_bs[0] = density->quadratic_bs[1];
    density->quadratic_cs[0] = density->quadratic_cs[1];

    density->quadratic_as[n - 1] = density->quadratic_as[n - 2];
    density->quadratic_bs[n - 1] = density->quadratic_bs[n - 2];
    density->quadratic_cs[n - 1] = density->quadratic_cs[n - 2];
}

void calculate_measurement_densities(beacon_localize_state_t *state)
{
    double *normalized_ys = NULL;

    for (int i = 0; i < state->n_tprops; i++) {
        all_measurements_t *all = &state->all_measurements[i];
        measurement_density_t *density = &state->densities[i];

        double spacing = PEAK_DELTA / 2;

        double min_x = DBL_MAX;
        double max_x = -DBL_MAX;

        double mean = 0;
        // double sum_x = 0;
        for (int j = 0; j < all->n; j++) {
            double v = all->tprops[j];
            if (v == 1337.42) {
                printf("Exit D\n");
                exit(1);
            }
            if (v < min_x) {
                min_x = v;
            }
            if (all->tprops[j] > max_x) {
                max_x = v;
            }
            mean += v;
        }
        mean /= all->n;

        double sigma = 1.3e-4;
        double inv_sigma = 1.0 / sigma;

        min_x -= sigma * 6;
        max_x += sigma * 6;

        int n = (int)((max_x - min_x) / spacing + 0.5) + 1;
        density->n = n;
        density->dx = spacing;
        density->min_x = min_x;
        density->ys = realloc(density->ys, sizeof(double) * n);
        density->mean_x = mean;

        // double a = (x - mu) / sigma;
        // return 1.0 / (sigma * sqrt(2 * M_PI)) * exp(-0.5 * a * a);

        double y_sum = 0;

        for (int j = 0; j < n; j++) {
            double x = min_x + spacing * j;
            double prob = eval_gauss_peak_all(all->n, x, all->tprops, inv_sigma);
            density->ys[j] = prob;
            y_sum += prob;
        }

        // normalized_ys is used to make the sampling table
        normalized_ys = realloc(normalized_ys, sizeof(double) * n);
        for (int j = 0; j < n; j++) {
            normalized_ys[j] = density->ys[j] / y_sum;
            // double x = min_x + spacing * j;
            // printf("%.8e %.8e\n", x, density->ys[j]);
        }
        // exit(1);

        // for (int j = 0; j < n; j++) {
        //     density->ys[j] = fastlog(density->ys[j]);
        // }

        calculate_qudratic_density_fit(density);

        // int zoom_n = n * 4;
        // double zoom_spacing = spacing * n / (double)zoom_n;
        // for (int j = 0; j < zoom_n; j++) {
        //     double x = min_x + zoom_spacing * j;
        //     double loc = (x - min_x) / density->dx;
        //     int best_i = (int)(loc + 0.5);
        //     double a = density->quadratic_as[best_i];
        //     double b = density->quadratic_bs[best_i];
        //     double c = density->quadratic_cs[best_i];
        //     double p = a * x * x + b * x + c;
        //     // double p2 = density->ys[best_i];
        //     printf("%.8e %.8e\n", x, p);
        // }
        // exit(1);

        // A length of 1000 will allow sampling the points down to 0.1% chance
        int len = 1000;
        density->sampling_len = len;
        density->sampling_table = malloc(len * sizeof(double));

        int table_i = 0;
        double prob_sum = 0;
        for (int j = 0; j < density->n; j++) {
            prob_sum += normalized_ys[j];
            while ((prob_sum >= (table_i + 1.0) / len) ||
                   (j == density->n - 1 && table_i < len)) {
                density->sampling_table[table_i] = density->min_x + j * density->dx;
                // printf("%d: %.5f\n", table_i, density->sampling_table[table_i]);
                table_i++;
            }
        }
    }

    free(normalized_ys);
}

void calc_stats(beacon_localize_state_t *state, all_measurements_t *all, double *variance, double *median)
{
    int n = all->n;

    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += all->tprops[i];
    }
    double avg = sum / n;

    sum = 0;
    for (int i = 0; i < n; i++) {
        double diff = all->tprops[i] - avg;
        sum += diff * diff;
    }
    *variance = sum / n;
    *median = sort_and_median(all);
}

// returns true if any outliers were pruned
bool prune_outliers(beacon_localize_state_t *state, all_measurements_t *all, double variance, double median)
{
    // remove points more than several standard deviations away from median
    // we want to get rid of true outliers
    double inv_stdev = 1.0 / sqrt(variance);
    int n_removed = 0;
    for (int i = 0; i < all->n; i++) {
        double val = all->tprops[i];
        double normalized = fabs(val - median) * inv_stdev;
        if (val < 0.5 || normalized > 4.0) {
            // swap remove
            all->tprops[i] = all->tprops[all->n - 1];
            all->n--;
            i--;
            n_removed++;
        }
    }
    state->measurements_removed += n_removed;
    return n_removed > 0;
}

void measurements_to_peaks(beacon_localize_state_t *state, all_measurements_t *all, zarray_t *peaks, double *overall_variance, int *modes_n, double *check)
{
    double variance = 0;
    double median = 0;
    calc_stats(state, all, &variance, &median);
    // printf("n: %d variance: %f median: %f\n", all->n, variance, median);

    int last_peaks_n = zarray_size(peaks);
    t_prop_peak_t last_data[last_peaks_n];
    memcpy(last_data, peaks->data, sizeof(last_data));

    find_peaks(all->n, all->tprops, variance, peaks, modes_n);
    int peaks_n = zarray_size(peaks);
    double min_diff = DBL_MAX;

    if (peaks_n == 0) {
        exit(1);
    }

    t_prop_peak_t* data = (t_prop_peak_t*)peaks->data;
    for (int i = 0; i < peaks_n; i++) {
        for (int j = 0; j < last_peaks_n; j++) {
            double diff = fabs(data[i].tprop - last_data[j].tprop);
            if (diff < min_diff) {
                min_diff = diff;
            }
        }
    }

    if (state->use_whole_measurement_variance) {
        for (int i = 0; i < peaks_n; i++) {
            data[i].variance = variance;
        }
    }

    if (overall_variance) {
        *overall_variance = variance;
    }
    if (check) {
        *check = min_diff;
    }
}

void process_measurements(beacon_localize_state_t *state, int measure_idx)
{
    all_measurements_t *all = &state->all_measurements[measure_idx];

    if (state->use_measurement_density) {
        double variance = 0;
        double median = 0;
        calc_stats(state, all, &variance, &median);
        bool continue_pruning = true;
        while (continue_pruning) {
            continue_pruning = prune_outliers(state, all, variance, median);
            if (continue_pruning) {
                calc_stats(state, all, &variance, &median);
            }
        }

        if (!state->silent_output) {
            int id_a = state->ids[state->as[measure_idx]];
            int id_b = state->ids[state->bs[measure_idx]];
            printf("%3d, %3d with %4d measurements\n", id_a, id_b, all->n);
        }
    } else {
        double variance = 0;
        double median = 0;
        calc_stats(state, all, &variance, &median);
        bool continue_pruning = true;
        while (continue_pruning) {
            continue_pruning = prune_outliers(state, all, variance, median);
            if (continue_pruning) {
                calc_stats(state, all, &variance, &median);
            }
        }

        if (!state->silent_output) {
            int modes_n = 0;
            double check = 0;
            zarray_t *peaks = state->all_t_prop_peaks[measure_idx];
            measurements_to_peaks(state, all, peaks, &variance, &modes_n, &check);
            t_prop_peak_t* data = (t_prop_peak_t*)peaks->data;
            int peaks_n = zarray_size(peaks);
            int id_a = state->ids[state->as[measure_idx]];
            int id_b = state->ids[state->bs[measure_idx]];
            printf("%3d, %3d reached variance %11.4e, change %11.4e, with %4d measurements, %d modes, %d peaks, all t_props ",
                   id_a, id_b, variance, check, all->n, modes_n, peaks_n);
            for (int i = 0; i < peaks_n; i++) {
                printf("%.7f (%.1e, %.2f), ", data[i].tprop, data[i].variance, data[i].weight);
            }
            printf("\n");
        }
    }
}

double eval_rmse_from_coords(beacon_localize_state_t *state, double *coords);

void calc_normal(double *coords, int *is, double *normal)
{
    // cross propuct between points
    double v1[3] = { coords[is[0] * 3 + 0] - coords[is[1] * 3 + 0],
                     coords[is[0] * 3 + 1] - coords[is[1] * 3 + 1],
                     coords[is[0] * 3 + 2] - coords[is[1] * 3 + 2] };

    double v2[3] = { coords[is[0] * 3 + 0] - coords[is[2] * 3 + 0],
                     coords[is[0] * 3 + 1] - coords[is[2] * 3 + 1],
                     coords[is[0] * 3 + 2] - coords[is[2] * 3 + 2] };

    normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
    normal[1] = -v1[0] * v2[2] + v1[2] * v2[0];
    normal[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void resolve_symmetry_ambiguity(beacon_localize_state_t *state, double *coords)
{
    // look for how many locations we have of the self_id, as from the "robot dance"
    int self_ids[3];
    int self_i = 0;
    bool skipped_second = false;
    for (int i = 0; i < state->n_locs && self_i < 3; i++) {
        int id = state->ids[i];
        if (id == state->self_id) {
            if (self_i == 1 && !skipped_second) {
                // skip second self beacon location which might be in a straight line with 0 and 2
                skipped_second = true;
            } else {
                self_ids[self_i] = i;
                self_i++;
            }
        }
    }

    if (self_i >= 3) {
        // We can use the calibration dance points to resolve the symmetry ambiguity. We know the dance
        // goes counter-clockwise, so we just need to make sure the cross-product also gives us that.
        double normal[3];
        double normal_coords[9] = { coords[self_ids[0] * 2], coords[self_ids[0] * 2 + 1], 0,
                                    coords[self_ids[1] * 2], coords[self_ids[1] * 2 + 1], 0,
                                    coords[self_ids[2] * 2], coords[self_ids[2] * 2 + 1], 0};
        int normal_is[3] = { 0, 1, 2 };
        calc_normal(normal_coords, normal_is, normal);
        bool flip_x = normal[2] < 0;
        for (int i = 0; i < state->n_locs; i++) {
            coords[i * 2] *= (flip_x ? -1 : 1);
        }
    }
}

double update_estimation(beacon_localize_state_t *state)
{
    if (state->use_measurement_density) {
        calculate_measurement_densities(state);
    }

    int n_locs = state->n_locs;
    int n_delays = state->n_beacons;
    int n_tprops = state->n_tprops;
    int n_coords = n_locs * 2;
    double coords[n_coords];

    int z_plane_delay_idx = state->ds[state->beacon_loc_i[state->self_id]];

    double residual = DBL_MAX;

    if (state->use_measurement_density) {
        calib_data_t *cal_data = beacon_autocal_create(n_locs, n_delays, n_tprops,
                                                       z_plane_delay_idx, state->robot_z_offset,
                                                       state->densities, state->as, state->bs, state->ds,
                                                       state->delay_mean, state->delay_sigma, state->recalibrate_delays);

        double best_delays[n_delays];
        double best_coords[n_locs * 2];

        for (int full_i = 0; full_i < state->n_full_rounds; full_i++) {
            double round_residual = DBL_MAX;
            double round_best_delays[n_delays];
            double round_best_coords[n_locs * 2];
            double round_best_initial_tprops[n_tprops];

            for (int round_i = 0; round_i < state->n_perturbation_rounds; round_i++) {
                double initial_tprops[n_tprops];
                for (int i = 0; i < n_tprops; i++) {
                    measurement_density_t *density = &state->densities[i];
                    // double rand_val = (double)rand() / RAND_MAX;
                    // double prob_sum = 0;
                    // for (int j = 0; j < density->n; j++) {
                    //     prob_sum += density->ys[j];
                    //     if (prob_sum >= rand_val) {
                    //         initial_tprops[i] = density->min_x + j * density->dx;
                    //         break;
                    //     }
                    // }
                    if (round_i == 0) {
                        initial_tprops[i] = density->mean_x;
                    } else {
                        int rand_i = rand() % density->sampling_len;
                        initial_tprops[i] = density->sampling_table[rand_i];
                    }
                    // double val2 = (double[]){0.530223357372, 0.525371844952, 0.530880659054, 0.518564077524, 0.530129457131, 0.521005483774, 0.535826071715, 0.521240234375, 0.535622621194, 0.535231370192, 0.525810046074, 0.520426432292, 0.527625450721, 0.525966546474, 0.526076096755, 0.536499023438, 0.522523537660, 0.531021509415, 0.532523913261, 0.520222981771, 0.526623848157, 0.526780348558, 0.525371844952, 0.522226186899, 0.524166791867, 0.529503455529, 0.524119841747, 0.525340544872}[i];
                    // initial_tprops[i] = val2;
                }

                double new_residual = beacon_autocal_trilaterate(cal_data, NULL, initial_tprops, state->delays, coords);

                if (fabs(new_residual) < 1e6) {
                    printf("\r%.6f", new_residual);
                }

                if (new_residual < round_residual) {
                    round_residual = new_residual;
                    // resolve_symmetry_ambiguity(state, coords);
                    // double rmse = eval_rmse_from_coords(state, coords);
                    // if (!state->silent_output) {
                    //     printf("\nrmse: %.6f\n", rmse);
                    // }

                    memcpy(round_best_delays, state->delays, sizeof(round_best_delays));
                    memcpy(round_best_coords, coords, sizeof(round_best_coords));
                    memcpy(round_best_initial_tprops, initial_tprops, sizeof(round_best_initial_tprops));
                }
            }

            if (!state->triangulation_only) {
                round_residual = beacon_autocalibrate(cal_data, NULL, round_best_initial_tprops, round_best_delays, round_best_coords, false);
            }

            if (round_residual < residual) {
                residual = round_residual;
                memcpy(best_delays, round_best_delays, sizeof(best_delays));
                memcpy(best_coords, round_best_coords, sizeof(best_coords));
            }
        }

        beacon_autocal_destroy(cal_data);

        if (!state->silent_output) {
            printf("\n      FINAL    %.6f\n", residual);
        }

        memcpy(state->delays, best_delays, sizeof(best_delays));
        memcpy(coords, best_coords, sizeof(best_coords));
    } else {
        zarray_t **all_peaks_to_use = state->all_t_prop_peaks;

        calib_data_t *cal_data = beacon_autocal_create(n_locs, n_delays, n_tprops,
                                                       z_plane_delay_idx, state->robot_z_offset,
                                                       NULL, state->as, state->bs, state->ds,
                                                       state->delay_mean, state->delay_sigma, state->recalibrate_delays);

        // derive all the peaks from the data
        for (int i = 0; i < n_tprops; i++) {
            all_measurements_t *all = &state->all_measurements[i];
            int modes_n = 0;
            double check = 0;
            double variance = 0;
            zarray_t *peaks = state->all_t_prop_peaks[i];
            measurements_to_peaks(state, all, peaks, &variance, &modes_n, &check);
            // int peaks_n = zarray_size(peaks);
            // t_prop_peak_t *peak_data = (t_prop_peak_t*)peaks->data;
            // printf("%3d, %3d has variance %11.4e with %4d measurements, %d modes, %d peaks, all t_props ",
            //        state->ids[state->as[i]], state->ids[state->bs[i]], variance, all->n, modes_n, peaks_n);
            // for (int i = 0; i < peaks_n; i++) {
            //     printf("%.7f (%.1e, %.2f), ", peak_data[i].tprop, peak_data[i].variance, peak_data[i].weight);
            // }
            // printf("\n");
        }

        if (state->try_alternate_tprops || state->first_peak || state->max_peak) {
            all_peaks_to_use = malloc(sizeof(zarray_t*) * n_tprops);
            for (int i = 0; i < n_tprops; i++) {
                zarray_t *peaks = state->all_t_prop_peaks[i];
                t_prop_peak_t *peak_data = (t_prop_peak_t*)peaks->data;
                all_peaks_to_use[i] = zarray_create(sizeof(t_prop_peak_t));

                if (state->max_peak) {
                    int peaks_n = zarray_size(peaks);
                    double max_weight = -DBL_MAX;
                    int max_idx = 0;
                    for (int j = 0; j < peaks_n; j++) {
                        if (peak_data[j].weight > max_weight) {
                            max_weight = peak_data[j].weight;
                            max_idx = j;
                        }
                    }
                    if (peak_data[max_idx].variance > 1e-6) {
                        peak_data[max_idx].variance = 1.7e-8; // SEE THE CRAZY IMPROVEMENT THIS MAKES!!!
                    }
                    zarray_add(all_peaks_to_use[i], &peak_data[max_idx]);
                    // double val2 = (double[]){0.5302234, 0.5253718, 0.5308807, 0.5185641, 0.5301295, 0.5210055, 0.5358261, 0.5212402, 0.5356226, 0.5352314, 0.5258100, 0.5204264, 0.5276255, 0.5259665, 0.5260761, 0.5364990, 0.5225235, 0.5310215, 0.5325239, 0.5202230, 0.5266238, 0.5267803, 0.5253718, 0.5222262, 0.5241668, 0.5295035, 0.5241198, 0.5253405}[i];
                    // t_prop_peak_t peak = {.tprop = val2, .variance = 1.7e-8, .weight = 0, .log_variance = -7.47};
                    // // peak.tprop = peak_data[max_idx].tprop;
                    // // peak.weight = peak_data[max_idx].weight;
                    // peak.variance = peak_data[max_idx].variance;
                    // // peak.log_variance = peak_data[max_idx].log_variance;
                    // zarray_add(all_peaks_to_use[i], &peak);
                } else {
                    zarray_add(all_peaks_to_use[i], &peak_data[0]);
                }
            }

            if (state->triangulation_only) {
                residual = beacon_autocal_trilaterate(cal_data, all_peaks_to_use, NULL, state->delays, coords);
            } else {
                residual = beacon_autocalibrate(cal_data, all_peaks_to_use, NULL, state->delays, coords, false);
            }

            if (!state->first_peak && !state->max_peak) {
                if (!state->silent_output) {
                    printf("\ninitial residual %.6f\n", residual);
                }

                int best_peaks[n_tprops];
                memset(best_peaks, 0, sizeof(best_peaks));

                for (int round_i = 0; round_i < 3; round_i++) {
                    for (int i = 0; i < n_tprops; i++) {
                        zarray_t *peaks = state->all_t_prop_peaks[i];
                        t_prop_peak_t *peak_data = (t_prop_peak_t*)peaks->data;
                        int peak_n = zarray_size(peaks);
                        if (peak_n <= 1) {
                            continue;
                        }

                        for (int j = 0; j < peak_n; j++) {
                            zarray_set(all_peaks_to_use[i], 0, &peak_data[j], NULL);
                            if (fabs(peak_data[j].tprop - peak_data[best_peaks[i]].tprop) <= ((double[]){ 1e-2, 1e-3, 1e-4, 1e-5})[round_i]) {
                                // if (fabs(peak_data[j].weight - peak_data[best_peaks[i]].weight) < 10) {
                                    continue;
                                // }
                            }

                            // printf("%e\n",    fabs(peak_data[j] - peak_data[0]));
                            double new_residual = beacon_autocalibrate(cal_data, all_peaks_to_use, NULL, state->delays, coords, false);
                            // if (!state->silent_output) {
                            //     printf("i %2d j %d new_residual %.6e\n", i, j, new_residual);
                            // }

                            if (new_residual < residual) {
                                residual = new_residual;
                                best_peaks[i] = j;
                            }
                        }
                        zarray_set(all_peaks_to_use[i], 0, &peak_data[best_peaks[i]], NULL);
                    }
                    if (!state->silent_output) {
                        printf("\nround %2d residual %.6f\n", round_i, residual);
                    }
                }
            }
        } else {
            beacon_autocalibrate(cal_data, all_peaks_to_use, NULL, state->delays, coords, false);
        }

        if (!state->silent_output) {
            beacon_autocalibrate(cal_data, all_peaks_to_use, NULL, state->delays, coords, true);
        }

        beacon_autocal_destroy(cal_data);

        if (all_peaks_to_use != state->all_t_prop_peaks) {
            for (int i = 0; i < n_tprops; i++) {
                zarray_destroy(all_peaks_to_use[i]);
            }
            free(all_peaks_to_use);
        }
    }

    for (int i = 0; i < n_locs; i++) {
        state->xs[i] = coords[i * 2];
        state->ys[i] = coords[i * 2 + 1];
    }

    if (!state->silent_output) {
        printf("\nDelays:\n");
        print_vec(n_delays, n_delays, state->delays);

        printf("Coords:\n");
        print_vec(n_coords, 2, coords);
    }

    state->estimation_ready = false;
    state->estimation_complete = true;

    return residual;
}

void clear_state(beacon_localize_state_t *state)
{
    free(state->as);
    free(state->bs);
    free(state->ds);
    free(state->ids);
    free(state->delays);
    free(state->xs);
    free(state->ys);

    for (int i = 0; i < state->n_tprops; i++) {
        zarray_destroy(state->all_t_prop_peaks[i]);
        free(state->densities[i].ys);
    }
    free(state->all_t_prop_peaks);
    free(state->all_measurements);
    free(state->densities);

    memset(state, 0, sizeof(*state));
}

void init_eval_los(beacon_localize_state_t *state);
void init_eval_nlos(beacon_localize_state_t *state);
void init_eval_nlos2(beacon_localize_state_t *state);
void init_eval_nlos3(beacon_localize_state_t *state);

void init_state(beacon_localize_state_t *state)
{
    double default_antenna_delay = DEFAULT_ANTENNA_DELAY;

    state->self_id = -1;
    state->recalibrate_delays = true;

    state->delay_mean = default_antenna_delay;
    state->delay_sigma = 3.3e-4;

    state->evaluate_by_full_rounds = getopt_get_bool(state->gopt, "by-full-rounds");
    state->evaluate_by_pert_rounds = getopt_get_bool(state->gopt, "by-pert-rounds");
    state->evaluate_by_measurements = getopt_get_bool(state->gopt, "by-measurements");
    state->evaluate_by_dancing = getopt_get_bool(state->gopt, "by-dancing");
    state->evaluate_by_delay_mean = getopt_get_bool(state->gopt, "by-delay-mean");
    state->evaluate_by_delay_sigma = getopt_get_bool(state->gopt, "by-delay-sigma");

    state->eval_run_i = getopt_get_int(state->gopt, "eval-run-i");
    state->fixed_delays = getopt_get_bool(state->gopt, "fixed-delays");

    if (getopt_get_bool(state->gopt, "first-peak")) {
        state->first_peak = true;
    } else if (getopt_get_bool(state->gopt, "max-peak")) {
        state->max_peak = true;
    }
    if (getopt_get_bool(state->gopt, "triangulation")) {
        state->triangulation_only = true;
    }
    if (!state->first_peak && !state->max_peak && !getopt_get_bool(state->gopt, "old")) {
        state->use_measurement_density = true;
    }

    if (getopt_get_bool(state->gopt, "evaluate-los")) {
        init_eval_los(state);
    } else if (getopt_get_bool(state->gopt, "evaluate-nlos")) {
        init_eval_nlos(state);
    } else if (getopt_get_bool(state->gopt, "evaluate-nlos2")) {
        init_eval_nlos2(state);
    } else if (getopt_get_bool(state->gopt, "evaluate-nlos3")) {
        init_eval_nlos3(state);
    }

    for (int i = 0; i < MAX_BEACONS; i++) {
        state->beacon_loc_i[i] = -1;
    }
}

typedef struct {
    int n;
    double *xs1;
    double *ys1;
    double *xs2;
    double *ys2;
} coords_match_data_t;

double coords_match_f(double *xyt, void *user)
{
    coords_match_data_t *data = user;
    int n = data->n;

    double xs1_prime[n];
    double ys1_prime[n];

    double cos_th = cos(xyt[2]);
    double sin_th = sin(xyt[2]);
    for (int i = 0; i < n; i++) {
        double x = data->xs1[i];
        double y = data->ys1[i];
        xs1_prime[i] = cos_th * x - sin_th * y + xyt[0];
        ys1_prime[i] = sin_th * x + cos_th * y + xyt[1];
    }

    double sum = 0;
    for (int i = 0; i < n; i++) {
        double dx = xs1_prime[i] - data->xs2[i];
        double dy = ys1_prime[i] - data->ys2[i];
        // exclude those points marked with NAN
        if (isfinite(dx) && isfinite(dy)) {
            sum += dx * dx + dy * dy;
        }
    }

    return sum;
}

double read_f64(FILE *f)
{
    double val;
    fread((void*)(&val), sizeof(val), 1, f);
    return val;
}

double read_f32(FILE *f)
{
    float val;
    fread((void*)(&val), sizeof(val), 1, f);
    return val;
}

int read_u8(FILE *f)
{
    uint8_t val;
    fread((void*)(&val), sizeof(val), 1, f);
    return val;
}

int read_u16(FILE *f)
{
    uint16_t val;
    fread((void*)(&val), sizeof(val), 1, f);
    return val;
}

// for debugging
void spoof_data_from_raw_dataset_file(beacon_localize_state_t *state, const char *file, bool use_dance_locations)
{
    FILE *f = fopen(file, "rb");
    if (!f) {
        fprintf(stderr, "Could not open %s\n", file);
        exit(1);
    }

    int n_beacons = 8;
    int *ids = (int[]){ 8, 12, 14, 16, 54, 61, 81, 87 };

    for (int i = 0; i < n_beacons; i++) {
        add_new_loc(state, ids[i]);
    }
    if (use_dance_locations) {
        add_new_loc(state, state->self_id);
        add_new_loc(state, state->self_id);
        add_new_loc(state, state->self_id);
        add_new_loc(state, state->self_id);
        add_new_loc(state, state->self_id);
    }

    for (int i = 0; i < state->n_beacons; i++) {
        state->delays[i] = state->delay_mean;
    }

    int measure_idx = 0;
    while (1) {
        double tprop_us = read_f64(f);
        int id_init = read_u8(f);
        int id_resp = read_u8(f);
        if (feof(f)) {
            break;
        }

        for (int i = 0; i < 3; i++) {
            read_f32(f); // temperature
            read_u16(f); // pp_amp
            read_u16(f); // impulse_pwr
            read_u16(f); // std_noise
            read_f32(f); // fp_power
            read_f32(f); // rx_power
        }

        for (int i = measure_idx; i < state->n_tprops; i++) {
            int a = state->as[i];
            int b = state->bs[i];
            int id_a = state->ids[a];
            int id_b = state->ids[b];
            if ((id_a == id_init && id_b == id_resp) ||
                (id_b == id_init && id_a == id_resp)) {
                measure_idx = i;
                break;
            }
        }

        all_measurements_t *all = &state->all_measurements[measure_idx];
        if (all->n < MAX_MEASUREMENTS) {
            int i = all->n;
            all->tprops[i] = tprop_us;
            all->n++;
        }
    }

    for (int i = 0; i < state->n_tprops; i++) {
        process_measurements(state, i);
    }

    fclose(f);
}

// from https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
double gaussian_random(double sigma)
{
    static const double epsilon = DBL_MIN;
    static const double two_pi = 2.0*3.14159265358979323846;

    static double z1;
    static bool generate = false;
    generate = !generate;
    if (!generate) {
        return z1 * sigma;
    }

    double u1, u2;
    do {
        u1 = rand() * (1.0 / RAND_MAX);
        u2 = rand() * (1.0 / RAND_MAX);
    } while (u1 <= epsilon);

    double z0;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma;
}

// with ground truth in feet
void init_eval(beacon_localize_state_t *state, int n_locs, double *ground_truth_xs, double *ground_truth_ys, bool ft_to_meters)
{
    double min_x = DBL_MAX;
    double max_x = -DBL_MAX;
    double min_y = DBL_MAX;
    double max_y = -DBL_MAX;

    // to meters
    for (int i = 0; i < n_locs; i++) {
        double x = ground_truth_xs[i];
        double y = ground_truth_ys[i];
        if (ft_to_meters) {
            x *= 12 * 0.0254;
            y *= 12 * 0.0254;
        }

        if (x < min_x) {
            min_x = x;
        }
        if (x > max_x) {
            max_x = x;
        }
        if (y < min_y) {
            min_y = y;
        }
        if (y > max_y) {
            max_y = y;
        }

        state->ground_truth_xs[i] = x;
        state->ground_truth_ys[i] = y;
    }

    state->eval_min_x = min_x;
    state->eval_max_x = max_x;
    state->eval_min_y = min_y;
    state->eval_max_y = max_y;
}

void init_eval_los(beacon_localize_state_t *state)
{
    // in feet
    double ground_truth_xs[] = { 0, 15, 2, NAN, 20, 10, 17, 8 };
    double ground_truth_ys[] = { 0,  1, 9, NAN,  3,  0,  7, 8 };

    init_eval(state, 8, ground_truth_xs, ground_truth_ys, true);
    state->self_id = 16;
    state->robot_z_offset = 0.3112;
    state->use_dance_locations = true;
    state->eval_dataset = 1; //11
    state->recalibrate_delays = !state->fixed_delays;
}

void init_eval_nlos(beacon_localize_state_t *state)
{
    // in feet, for dataset 12
    // double ground_truth_xs[] = { 24.333, 3,     10, NAN,      -10,  4, -10, 28.167 };
    // double ground_truth_ys[] = { 17,     8, -12.25, NAN, -12.9167, 24,  25,    -13 };

    // for dataset 14 (same as los positions)
    double ground_truth_xs[] = { 0, 15, 2, NAN, 20, 10, 17, 8 };
    double ground_truth_ys[] = { 0,  1, 9, NAN,  3,  0,  7, 8 };

    init_eval(state, 8, ground_truth_xs, ground_truth_ys, true);
    state->self_id = 16;
    state->robot_z_offset = 0.3112;
    state->use_dance_locations = true;
    state->eval_dataset = 2;
    state->recalibrate_delays = !state->fixed_delays;
}

void init_eval_nlos2(beacon_localize_state_t *state)
{
    // in feet
    double ground_truth_xs[] = { 0, 15, 2, NAN, 20, 10, 17, 8 };
    double ground_truth_ys[] = { 0,  1, 9, NAN,  3,  0,  7, 8 };

    init_eval(state, 8, ground_truth_xs, ground_truth_ys, true);
    state->self_id = 16;
    state->robot_z_offset = 0.3112;
    state->use_dance_locations = true;
    state->eval_dataset = 3;
    state->recalibrate_delays = !state->fixed_delays;
}

void init_eval_nlos3(beacon_localize_state_t *state)
{
    // in feet
    double ground_truth_xs[] = { 0, 15, 2, NAN, 20, 10, 17, 8 };
    double ground_truth_ys[] = { 0,  1, 9, NAN,  3,  0,  7, 8 };

    init_eval(state, 8, ground_truth_xs, ground_truth_ys, true);
    state->self_id = 16;
    state->robot_z_offset = 0.3112;
    state->use_dance_locations = false;
    state->eval_dataset = 4;
    state->recalibrate_delays = !state->fixed_delays;
}

double inner_eval_rmse_xyt(beacon_localize_state_t *state, double *xs, double *ys, double *xyt)
{
    // printf("\n%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", xs[0], xs[1], xs[2], xs[3], xs[4], xs[5], xs[6], xs[7]);

    int eval_locs = 8;

    coords_match_data_t data = { .n = eval_locs, .xs1 = xs, .ys1 = ys,
                                 .xs2 = state->ground_truth_xs, .ys2 = state->ground_truth_ys };

    // Now find the best rotation between the two
    xyt[0] = 0;
    xyt[1] = 0;
    xyt[2] = 0;
    double best_residual = DBL_MAX;
    double best_xyt[3] = { 0, 0, 0 };
    for (int i = -18; i < 18; i++) {
        xyt[2] = i * M_PI / 36.0;
        double residual = coords_match_f(xyt, &data);
        if (residual < best_residual) {
            best_residual = residual;
            best_xyt[2] = xyt[2];
        }
    }
    xyt[2] = best_xyt[2];

    // And then also the best x-y translation
    best_residual = newton_optimize(3, coords_match_f, &data, xyt);
    double rmse = sqrt(best_residual / eval_locs);
    return rmse;
}

double eval_rmse_xyt(beacon_localize_state_t *state, double *xs, double *ys, double *xyt)
{
    int eval_locs = 8;

    double res1 = inner_eval_rmse_xyt(state, xs, ys, xyt);

    // also try the mirror, if evaluating dancing and can't correctly resolve the symmetry
    // if (state->evaluate_by_dancing) {
        double mirror_ys[eval_locs];
        for (int i = 0; i < eval_locs; i++) {
            mirror_ys[i] = -ys[i];
        }
        double xyt2[3];
        double res2 = inner_eval_rmse_xyt(state, xs, mirror_ys, xyt2);
        if (res2 < res1) {
            for (int i = 0; i < eval_locs; i++) {
                ys[i] = mirror_ys[i];
            }
            memcpy(xyt, xyt2, sizeof(xyt2));
            return res2;
        }
    // }
    return res1;
}

double eval_rmse_from_coords(beacon_localize_state_t *state, double *coords)
{
    int eval_locs = 8;

    double xs[eval_locs];
    double ys[eval_locs];
    for (int i = 0; i < eval_locs; i++) {
        xs[i] = coords[i * 2];
        ys[i] = coords[i * 2 + 1];
    }

    double xyt[3];
    return eval_rmse_xyt(state, xs, ys, xyt);
}

void init_evaluate(beacon_localize_state_t *state)
{
    getopt_t *gopt = state->gopt;
    clear_state(state);
    state->gopt = gopt;

    init_state(state);

    state->try_alternate_tprops = true;
    state->silent_output = true;
    bool use_dance_locations = state->use_dance_locations;

    char file[256];
    snprintf(file, sizeof(file), "datasets/dataset%d.txt.binary", state->eval_dataset);
    spoof_data_from_raw_dataset_file(state, file, use_dance_locations);

    state->estimation_ready = true;
}

double inner_evaluate(beacon_localize_state_t *state, double *residual)
{
    *residual = update_estimation(state);
    // printf("\n");

    double xyt[3];
    double rmse = eval_rmse_xyt(state, state->xs, state->ys, xyt);
    // printf("xyt: %6.3f, %6.3f, %6.2f\t\tRMSE: %8.4f\n", xyt[0], xyt[1], xyt[2], rmse);

    double cos_th = cos(xyt[2]);
    double sin_th = sin(xyt[2]);
    for (int i = 0; i < state->n_locs; i++) {
        double x = state->xs[i];
        double y = state->ys[i];
        state->xs[i] = cos_th * x - sin_th * y + xyt[0];
        state->ys[i] = sin_th * x + cos_th * y + xyt[1];
        // state->xs[i] = state->ground_truth_xs[i];
        // state->ys[i] = state->ground_truth_ys[i];
    }

    // printf("\nCoords:\n");
    // for (int i = 0; i < state->n_locs; i++) {
    //     printf("%12.2f %12.2f\n", state->xs[i], state->ys[i]);
    // }

    return rmse;
}

void do_evaluate(beacon_localize_state_t *state)
{
    init_evaluate(state);

    int full_rounds_start = 20;
    int full_rounds_end = 21;
    if (state->evaluate_by_full_rounds) {
        full_rounds_start = 1;
        full_rounds_end = 11;
    }

    // overriding this for a clean 250 rounds below
    int pert_rounds_start = 9; //9;
    int pert_rounds_end = 10; //10;
    if (state->evaluate_by_pert_rounds) {
        pert_rounds_start = 0;
        pert_rounds_end = 14;
    }

    int measurements_start = 256;
    int measurements_end = 256;
    if (state->evaluate_by_measurements) {
        measurements_start = 1;
        measurements_end = 1024;
    }

    int dancing_tprop_n[] = { 28, 35, 42, 49, 56, 63 };
    int dancing_start = 0;
    int dancing_end = 1;
    int dancing_total = 6;
    if (state->evaluate_by_dancing) {
        dancing_start = 0;
        dancing_end = dancing_total;
    }

    double delay_means[] = { 0, 0.51, 0.511, 0.512, 0.513, 0.514, 0.515, 0.516, 0.517, 0.518, 0.519, 0.52, 0.521, 0.522 };
    int mean_start = 7;
    int mean_end = 8;
    if (state->evaluate_by_delay_mean) {
        mean_start = 1;
        mean_end = 14;
    }

    double delay_sigmas[] = { 3.3e-4*128, 3.3e-4*64, 3.3e-4*32, 3.3e-4*16, 3.3e-4*8, 3.3e-4*4, 3.3e-4*2, 3.3e-4, 3.3e-4/2, 3.3e-4/4 };
    int sigma_start = 7;
    int sigma_end = 8;
    if (state->evaluate_by_delay_sigma) {
        sigma_start = 0;
        sigma_end = 10;
        // if (!state->evaluate_by_delay_mean) {
        //     mean_start = 0;
        //     mean_end = 1;
        // }
    }

    // clean out empty empty measurements
    printf("\n");
    for (int i = 0; i < state->n_tprops; i++) {
        if (state->all_measurements[i].n < 50) {
            printf("Removed sparse (%d)  measurements between %d, %d\n", state->all_measurements[i].n, state->ids[state->as[i]], state->ids[state->bs[i]]);

            state->n_tprops--;

            int last = state->n_tprops;
            state->as[i] = state->as[last];
            state->bs[i] = state->bs[last];

            all_measurements_t tmp;
            memcpy(&tmp, &state->all_measurements[i], sizeof(tmp));
            memcpy(&state->all_measurements[i], &state->all_measurements[last], sizeof(tmp));
            memcpy(&state->all_measurements[last], &tmp, sizeof(tmp));

            measurement_density_t dense;
            memcpy(&dense, &state->densities[i], sizeof(dense));
            memcpy(&state->densities[i], &state->densities[last], sizeof(dense));
            memcpy(&state->densities[last], &dense, sizeof(dense));

            zarray_t* peaks = state->all_t_prop_peaks[i];
            state->all_t_prop_peaks[i] = state->all_t_prop_peaks[last];
            state->all_t_prop_peaks[last] = peaks;

            for (int j = 0; j < dancing_end; j++) {
                if (i < dancing_tprop_n[j]) {
                    dancing_tprop_n[j]--;
                }
            }

            i--;
        }
    }

    int original_all_ns[state->n_tprops];
    for (int l = 0; l < state->n_tprops; l++) {
        original_all_ns[l] = state->all_measurements[l].n;
    }

    int original_n_tprops = state->n_tprops;
    int base_n_locs = state->n_locs - (dancing_total - 1);

    for (int sigma_i = sigma_start; sigma_i < sigma_end; sigma_i++) {
        state->delay_sigma = delay_sigmas[sigma_i];
        for (int mean_i = mean_start; mean_i < mean_end; mean_i++) {
            state->delay_mean = delay_means[mean_i];
            for (int dancing_i = dancing_start; dancing_i < dancing_end; dancing_i++) {
                if (dancing_tprop_n[dancing_i] < original_n_tprops) {
                    state->n_tprops = dancing_tprop_n[dancing_i];
                    state->n_locs = base_n_locs + dancing_i;
                }

                // for (int i = 0; i < 14; i++) {
                // for (int i = 6; i <= 6; i++) {
                for (int pert_i = pert_rounds_start; pert_i < pert_rounds_end; pert_i++) {
                    int n_perturbation_rounds = (int)(pow(10, pert_i * 0.25) + 0.5);
                    if (!state->evaluate_by_pert_rounds) {
                        n_perturbation_rounds = 250;
                    }
                    state->n_perturbation_rounds = n_perturbation_rounds;

                    for (int full_n = full_rounds_start; full_n < full_rounds_end; full_n++) {
                        state->n_full_rounds = full_n;

                        for (int measure_max = measurements_start; measure_max <= measurements_end; measure_max *= 2) {
                            for (int j = 0; j < 10; j++) {
                                srand(j + 1000 * state->eval_run_i);

                                double residual;
                                // randomly choose the samples to use for each link/edge
                                for (int l = 0; l < state->n_tprops; l++) {
                                    all_measurements_t *all = &state->all_measurements[l];
                                    int measurements_n = original_all_ns[l];
                                    if (measurements_n > measure_max) {
                                        measurements_n = measure_max;
                                    }
                                    all->n = measurements_n;
                                    for (int m = 0; m < measurements_n; m++) {
                                        int rand_i = rand() % original_all_ns[l];
                                        double val = all->tprops[m];
                                        all->tprops[m] = all->tprops[rand_i];
                                        all->tprops[rand_i] = val;
                                        if (all->tprops[m] == 1337.42) {
                                            printf("Exit C\n");
                                            exit(1);
                                        }
                                    }
                                }

                                uint64_t start_us = utime_now();
                                double rmse = inner_evaluate(state, &residual);
                                double timed_ms = (utime_now() - start_us) / 1000.0;
                                printf("With ");
                                if (state->evaluate_by_full_rounds) {
                                    printf("Full Rounds = %4d, ", full_n);
                                }
                                if (state->evaluate_by_pert_rounds) {
                                    printf("Pert Rounds = %4d, ", n_perturbation_rounds);
                                }
                                if (state->evaluate_by_measurements) {
                                    printf("Measurements = %4d, ", measure_max);
                                }
                                if (state->evaluate_by_dancing) {
                                    printf("Dancing = %d, ", dancing_i);
                                }
                                if (state->evaluate_by_delay_mean) {
                                    printf("Delay-Mean = %d, ", mean_i);
                                }
                                if (state->evaluate_by_delay_sigma) {
                                    printf("Delay-Sigma = %d, ", sigma_i);
                                }
                                printf("ms: %6.2f, res: %8.4f, %10.6f", timed_ms, residual, rmse);

                                if (state->evaluate_by_delay_mean) {
                                    for (int i = 0; i < state->n_beacons; i++) {
                                        printf(" %.8f", state->delays[i]);
                                    }
                                }

                                printf("\n");
                            }
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char **argv)
{
    beacon_localize_state_t state = { 0 };

    getopt_t *gopt = getopt_create();
    getopt_add_bool(gopt, 'h', "help", 0, "Show usage");
    getopt_add_bool(gopt, '\0', "debug", 0, "Debugging mode");
    getopt_add_bool(gopt, '\0', "evaluate-los", 0, "Do evaluation on los dataset");
    getopt_add_bool(gopt, '\0', "evaluate-nlos", 0, "Do evaluation on nlos dataset");
    getopt_add_bool(gopt, '\0', "evaluate-nlos2", 0, "Do evaluation on nlos2 dataset");
    getopt_add_bool(gopt, '\0', "evaluate-nlos3", 0, "Do evaluation on nlos3 dataset");
    getopt_add_bool(gopt, '\0', "by-full-rounds", 0, "Parameterize evaluation by number of full (including all perturbations) rounds");
    getopt_add_bool(gopt, '\0', "by-pert-rounds", 0, "Parameterize evaluation by number of random initialization (perturbation) rounds");
    getopt_add_bool(gopt, '\0', "by-measurements", 0, "Parameterize evaluation by number of measurements used");
    getopt_add_bool(gopt, '\0', "by-dancing", 0, "Parameterize evaluation by number of robot dance positions used");
    getopt_add_bool(gopt, '\0', "by-delay-mean", 0, "Parameterize evaluation by antenna delay prior mean");
    getopt_add_bool(gopt, '\0', "by-delay-sigma", 0, "Parameterize evaluation by antenna delay prior standard deviation");
    getopt_add_bool(gopt, '\0', "fixed-delays", 0, "Hold antenna delays to the nominal value; do not optimize");
    getopt_add_bool(gopt, '\0', "old", 0, "Do evaluation by old peak estimation method");
    getopt_add_bool(gopt, '\0', "first-peak", 0, "Do evaluation by first-peak method");
    getopt_add_bool(gopt, '\0', "max-peak", 0, "Do evaluation by max-peak method");
    getopt_add_bool(gopt, '\0', "triangulation", 0, "Do evaluation by triangulation-only method");
    getopt_add_int(gopt, '\0', "eval-run-i", "0", "Evaluation run iteration: evaluation results are deterministic for a give value of this");
    if (!getopt_parse(gopt, argc, argv, true) || getopt_get_bool(gopt, "help")) {
        getopt_do_usage(gopt);
        exit(getopt_get_bool(gopt, "help") ? EXIT_SUCCESS : EXIT_FAILURE);
    }
    state.gopt = gopt;

    if (getopt_get_bool(gopt, "evaluate-los") ||
        getopt_get_bool(gopt, "evaluate-nlos") ||
        getopt_get_bool(gopt, "evaluate-nlos2") ||
        getopt_get_bool(gopt, "evaluate-nlos3")) {

        do_evaluate(&state);
    }

    getopt_destroy(state.gopt);
    clear_state(&state);
    return 0;
}
