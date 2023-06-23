#include <math.h>

#include "summary_functions.h"

float
D(tsk_size_t w_AB, tsk_size_t w_Ab, tsk_size_t w_aB, tsk_size_t n)
{
    float p_AB = (float) w_AB / (float) n;
    float p_Ab = (float) w_Ab / (float) n;
    float p_aB = (float) w_aB / (float) n;

    float p_A = p_AB + p_Ab;
    float p_B = p_AB + p_aB;

    return p_AB - (p_A * p_B);
}

float
D2(tsk_size_t w_AB, tsk_size_t w_Ab, tsk_size_t w_aB, tsk_size_t n)
{
    return powf(D(w_AB, w_Ab, w_aB, n), 2.0);
}

float
r2(tsk_size_t w_AB, tsk_size_t w_Ab, tsk_size_t w_aB, tsk_size_t n)
{
    float p_AB = (float) w_AB / (float) n;
    float p_Ab = (float) w_Ab / (float) n;
    float p_aB = (float) w_aB / (float) n;

    float p_A = p_AB + p_Ab;
    float p_B = p_AB + p_aB;

    float D_ = p_AB - (p_A * p_B);
    float denom = p_A * p_B * (1 - p_A) * (1 - p_B);

    if (denom == 0 && D_ == 0) {
        return NAN; // TODO: what sort of value should I actually return here
    }

    return (D_ * D_) / denom;
}

float
D_prime(tsk_size_t w_AB, tsk_size_t w_Ab, tsk_size_t w_aB, tsk_size_t n)
{
    float p_AB = (float) w_AB / (float) n;
    float p_Ab = (float) w_Ab / (float) n;
    float p_aB = (float) w_aB / (float) n;

    float p_A = p_AB + p_Ab;
    float p_B = p_AB + p_aB;

    float D_ = p_AB - (p_A * p_B);
    if (D_ >= 0) {
        return D_ / TSK_MIN(p_A * (1 - p_B), (1 - p_A) * p_B);
    }
    return D_ / TSK_MIN(p_A * p_B, (1 - p_A) * (1 - p_B));
}

float
r(tsk_size_t w_AB, tsk_size_t w_Ab, tsk_size_t w_aB, tsk_size_t n)
{
    float p_AB = (float) w_AB / (float) n;
    float p_Ab = (float) w_Ab / (float) n;
    float p_aB = (float) w_aB / (float) n;

    float p_A = p_AB + p_Ab;
    float p_B = p_AB + p_aB;

    float D_ = p_AB - (p_A * p_B);
    float denom = p_A * p_B * (1 - p_A) * (1 - p_B);

    if (denom == 0 && D_ == 0) {
        return NAN; // TODO: what sort of value should I actually return here
    }

    return D_ / sqrtf(denom);
}

float
Dz(tsk_size_t w_AB, tsk_size_t w_Ab, tsk_size_t w_aB, tsk_size_t n)
{
    float p_AB = (float) w_AB / (float) n;
    float p_Ab = (float) w_Ab / (float) n;
    float p_aB = (float) w_aB / (float) n;

    float p_A = p_AB + p_Ab;
    float p_B = p_AB + p_aB;

    float D_ = p_AB - (p_A * p_B);

    return D_ * (1 - 2 * p_A) * (1 - 2 * p_B);
}

float
pi2(tsk_size_t w_AB, tsk_size_t w_Ab, tsk_size_t w_aB, tsk_size_t n)
{
    float p_AB = (float) w_AB / (float) n;
    float p_Ab = (float) w_Ab / (float) n;
    float p_aB = (float) w_aB / (float) n;

    float p_A = p_AB + p_Ab;
    float p_B = p_AB + p_aB;

    return p_A * (1 - p_A) * p_B * (1 - p_B);
}
