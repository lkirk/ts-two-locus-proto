import numpy as np


def D(w_AB, w_Ab, w_aB, n):
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    # p_ab = 1 - p_AB - p_Ab - p_aB

    # return p_AB * p_ab - p_Ab * p_aB
    return p_AB - (p_A * p_B)


def D2(w_AB, w_Ab, w_aB, n):
    return D(w_AB, w_Ab, w_aB, n) ** 2


def r2(w_AB, w_Ab, w_aB, n):
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    D_ = p_AB - (p_A * p_B)
    denom = p_A * p_B * (1 - p_A) * (1 - p_B)

    if denom == 0 and D_ == 0:
        return np.nan

    return (D_ * D_) / denom


def D_prime(w_AB, w_Ab, w_aB, n):
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    D_ = p_AB - (p_A * p_B)
    if D_ >= 0:
        return D_ / min(p_A * (1 - p_B), (1 - p_A) * p_B)
    return D_ / min(p_A * p_B, (1 - p_A) * (1 - p_B))


def r(w_AB, w_Ab, w_aB, n):
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    D_ = p_AB - (p_A * p_B)
    denom = p_A * p_B * (1 - p_A) * (1 - p_B)

    if denom == 0 and D_ == 0:
        return np.nan

    return D_ / np.sqrt(denom)


def Dz(w_AB, w_Ab, w_aB, n):
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    D_ = p_AB - (p_A * p_B)

    return D_ * (1 - 2 * p_A) * (1 - 2 * p_B)


def pi2(w_AB, w_Ab, w_aB, n):
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    return p_A * (1 - p_A) * p_B * (1 - p_B)
