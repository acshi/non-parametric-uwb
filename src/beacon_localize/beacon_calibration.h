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

#pragma once

#include <stdint.h>
#include <stdbool.h>
#include "common/zarray.h"

#include "shared_structs.h"

struct calib_data;
typedef struct calib_data calib_data_t;

double range_bias_at(double m);
double dist_between_simple(double *coords, int a, int b);
double dist_from_tof(double flight_time);
void print_vec(int n, int per_row, double *vec);
bool preconditioned_conjugate_gradient(int n, double *A, double *b, double epsilon, int kmax, double *x);
double newton_optimize(int n, double (*f)(double*, void*), void* user, double *x);
double newton_optimize_base(int n, double (*f)(double*, void*), void* user, double *x, double epsilon, int kmax);
calib_data_t *beacon_autocal_create(int n_locs, int n_delays, int n_tprops,
                                    int z_plane_delay_idx, double z_plane_2d_offset,
                                    measurement_density_t *densities, int *as, int *bs, int *ds,
                                    double delay_mean, double delay_sigma, bool recalibrate_delays);
void beacon_autocal_destroy(calib_data_t *data);
double beacon_autocal_trilaterate(calib_data_t *data, zarray_t **tprop_peaks, double *initial_tprops, double *delays, double *coords);
double beacon_autocalibrate(calib_data_t *data, zarray_t **tprop_peaks, double *initial_tprops, double *delays, double *coords, bool verbose);

double fastexp(double p);
double fastlog(double x);

void find_fastlog2_coefs();
void find_fastpow2_coefs();
