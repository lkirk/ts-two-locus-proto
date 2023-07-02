#pragma once
#include "prototype.h"

// copy this over from tskit to test integration
typedef struct {
    const tsk_id_t *sample_sets;
    tsk_size_t num_sample_sets;
    const tsk_size_t *sample_set_sizes;
    const tsk_id_t *set_indexes;
} sample_count_stat_params_t;

int D(tsk_size_t state_dim, const double *state, tsk_size_t result_dim, double *result,
    void *params);
int D2(tsk_size_t state_dim, const double *state, tsk_size_t result_dim, double *result,
    void *params);
int r2(tsk_size_t state_dim, const double *state, tsk_size_t result_dim, double *result,
    void *params);
int D_prime(tsk_size_t state_dim, const double *state, tsk_size_t result_dim,
    double *result, void *params);
int r(tsk_size_t state_dim, const double *state, tsk_size_t result_dim, double *result,
    void *params);
int Dz(tsk_size_t state_dim, const double *state, tsk_size_t result_dim, double *result,
    void *params);
int pi2(tsk_size_t state_dim, const double *state, tsk_size_t result_dim, double *result,
    void *params);
