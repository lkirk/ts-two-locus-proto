#include <math.h>

#include "summary_functions.h"

/* D(tsk_size_t w_AB, tsk_size_t w_Ab, tsk_size_t w_aB) */
int
D(tsk_size_t state_dim, const double *state, tsk_size_t TSK_UNUSED(result_dim),
    double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    double n;
    const double *state_row;
    for (tsk_size_t j = 0; j < state_dim; j++) {
        n = (double) args.sample_set_sizes[j];
        state_row = GET_2D_ROW(state, 3, j);
        double p_AB = state_row[0] / n;
        double p_Ab = state_row[1] / n;
        double p_aB = state_row[2] / n;

        double p_A = p_AB + p_Ab;
        double p_B = p_AB + p_aB;
        result[j] = p_AB - (p_A * p_B);
    }

    return 0;
}

int
D2(tsk_size_t state_dim, const double *state, tsk_size_t TSK_UNUSED(result_dim),
    double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    double n;
    const double *state_row;
    for (tsk_size_t j = 0; j < state_dim; j++) {
        n = (double) args.sample_set_sizes[j];
        state_row = GET_2D_ROW(state, 3, j);
        double p_AB = state_row[0] / n;
        double p_Ab = state_row[1] / n;
        double p_aB = state_row[2] / n;

        double p_A = p_AB + p_Ab;
        double p_B = p_AB + p_aB;
        result[j] = p_AB - (p_A * p_B);
        result[j] *= result[j];
    }

    return 0;
}

int
r2(tsk_size_t state_dim, const double *state, tsk_size_t TSK_UNUSED(result_dim),
    double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    double n;
    const double *state_row;
    for (tsk_size_t j = 0; j < state_dim; j++) {
        n = (double) args.sample_set_sizes[j];
        state_row = GET_2D_ROW(state, 3, j);
        double p_AB = state_row[0] / n;
        double p_Ab = state_row[1] / n;
        double p_aB = state_row[2] / n;

        double p_A = p_AB + p_Ab;
        double p_B = p_AB + p_aB;

        double D_ = p_AB - (p_A * p_B);
        double denom = p_A * p_B * (1 - p_A) * (1 - p_B);

        if (denom == 0 && D_ == 0) {
            result[j] = NAN; // TODO: what sort of value should I actually return here
        } else {
            result[j] = (D_ * D_) / denom;
        }
    }
    return 0;
}

int
D_prime(tsk_size_t state_dim, const double *state, tsk_size_t TSK_UNUSED(result_dim),
    double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    double n;
    const double *state_row;
    for (tsk_size_t j = 0; j < state_dim; j++) {
        n = (double) args.sample_set_sizes[j];
        state_row = GET_2D_ROW(state, 3, j);
        double p_AB = state_row[0] / n;
        double p_Ab = state_row[1] / n;
        double p_aB = state_row[2] / n;

        double p_A = p_AB + p_Ab;
        double p_B = p_AB + p_aB;

        double D_ = p_AB - (p_A * p_B);
        if (D_ >= 0) {
            result[j] = D_ / TSK_MIN(p_A * (1 - p_B), (1 - p_A) * p_B);
        } else {
            result[j] = D_ / TSK_MIN(p_A * p_B, (1 - p_A) * (1 - p_B));
        }
    }
    return 0;
}

int
r(tsk_size_t state_dim, const double *state, tsk_size_t TSK_UNUSED(result_dim),
    double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    double n;
    const double *state_row;
    for (tsk_size_t j = 0; j < state_dim; j++) {
        n = (double) args.sample_set_sizes[j];
        state_row = GET_2D_ROW(state, 3, j);
        double p_AB = state_row[0] / n;
        double p_Ab = state_row[1] / n;
        double p_aB = state_row[2] / n;

        double p_A = p_AB + p_Ab;
        double p_B = p_AB + p_aB;

        double D_ = p_AB - (p_A * p_B);
        double denom = p_A * p_B * (1 - p_A) * (1 - p_B);

        if (denom == 0 && D_ == 0) {
            result[j] = NAN; // TODO: what sort of value should I actually return here
        } else {
            result[j] = D_ / sqrt(denom);
        }
    }
    return 0;
}

int
Dz(tsk_size_t state_dim, const double *state, tsk_size_t TSK_UNUSED(result_dim),
    double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    double n;
    const double *state_row;
    for (tsk_size_t j = 0; j < state_dim; j++) {
        n = (double) args.sample_set_sizes[j];
        state_row = GET_2D_ROW(state, 3, j);
        double p_AB = state_row[0] / n;
        double p_Ab = state_row[1] / n;
        double p_aB = state_row[2] / n;

        double p_A = p_AB + p_Ab;
        double p_B = p_AB + p_aB;

        double D_ = p_AB - (p_A * p_B);

        result[j] = D_ * (1 - 2 * p_A) * (1 - 2 * p_B);
    }
    return 0;
}

int
pi2(tsk_size_t state_dim, const double *state, tsk_size_t TSK_UNUSED(result_dim),
    double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    double n;
    const double *state_row;
    for (tsk_size_t j = 0; j < state_dim; j++) {
        n = (double) args.sample_set_sizes[j];
        state_row = GET_2D_ROW(state, 3, j);
        double p_AB = state_row[0] / n;
        double p_Ab = state_row[1] / n;
        double p_aB = state_row[2] / n;

        double p_A = p_AB + p_Ab;
        double p_B = p_AB + p_aB;
        result[j] = p_A * (1 - p_A) * p_B * (1 - p_B);
    }
    return 0;
}
