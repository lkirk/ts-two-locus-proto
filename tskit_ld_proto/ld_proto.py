"""
Prototype functionality for computing two-locus statistics in tskit

To test the functionality, here's some examples:

>>> from tskit_ld_proto.summary_functions import D, r2, r, D_prime, Dz, pi2

>>> correlated = (np.array([[0, 1, 1, 0, 2, 2, 1, 0, 1], [1, 2, 2, 1, 0, 0, 2, 1, 2]]), np.array([3, 3]))
>>> correlated_biallelic = (np.array([[0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1]]), np.array([2, 2]))
>>> uncorrelated = (np.array([[0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]]), np.array([3, 3]))
>>> uncorrelated_biallelic = (np.array([[0, 0, 0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, 0, 1, 1]]), np.array([2, 2]))
>>> repulsion_biallelic = (np.array([[0, 0, 0, 0, 1, 1, 1, 1], [1, 1, 1, 1, 0, 0, 0, 0]]), np.array([2, 2]))

>>> _two_site_general_stat(*correlated, D, total, polarized=True)
[[ 0.05555556 -0.01851852]
 [-0.01851852  0.04320988]]
>>> _two_site_general_stat(*correlated_biallelic, D, total, polarized=True)
[[0.25 0.25]
 [0.25 0.25]]
>>> _two_site_general_stat(*uncorrelated, D, total, polarized=True)
[[0.05555556 0.        ]
 [0.         0.05555556]]
>>> _two_site_general_stat(*uncorrelated_biallelic, D, total, polarized=True)
[[0.25 0.  ]
 [0.   0.25]]
>>> _two_site_general_stat(*repulsion_biallelic, D, total, polarized=True)
[[ 0.25 -0.25]
 [-0.25  0.25]]

>>> _two_site_general_stat(*correlated, D, total, polarized=False)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*correlated_biallelic, D, total, polarized=False)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*uncorrelated, D, total, polarized=False)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*uncorrelated_biallelic, D, total, polarized=False)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*repulsion_biallelic, D, total, polarized=False)
[[0. 0.]
 [0. 0.]]

>>> _two_site_general_stat(*correlated, D_prime, hap_weighted, polarized=True)
[[0.66666667 0.44444444]
 [0.44444444 0.77777778]]
>>> _two_site_general_stat(*correlated_biallelic, D_prime, hap_weighted, polarized=True)
[[0.5 0.5]
 [0.5 0.5]]
>>> _two_site_general_stat(*uncorrelated, D_prime, hap_weighted, polarized=True)
[[0.66666667 0.        ]
 [0.         0.66666667]]
>>> _two_site_general_stat(*uncorrelated_biallelic, D_prime, hap_weighted, polarized=True)
[[0.5 0. ]
 [0.  0.5]]
>>> _two_site_general_stat(*repulsion_biallelic, D_prime, hap_weighted, polarized=True)
[[0.5 0. ]
 [0.  0.5]]

>>> _two_site_general_stat(*correlated, D_prime, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]
>>> _two_site_general_stat(*correlated_biallelic, D_prime, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]
>>> _two_site_general_stat(*uncorrelated, D_prime, hap_weighted, polarized=False)
[[1. 0.]
 [0. 1.]]
>>> _two_site_general_stat(*uncorrelated_biallelic, D_prime, hap_weighted, polarized=False)
[[1. 0.]
 [0. 1.]]
>>> _two_site_general_stat(*repulsion_biallelic, D_prime, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]

>>> _two_site_general_stat(*correlated, D2, total, polarized=True)
[[0.02758726 0.02453894]
 [0.02453894 0.03856119]]
>>> _two_site_general_stat(*correlated_biallelic, D2, total, polarized=True)
[[0.0625 0.0625]
 [0.0625 0.0625]]
>>> _two_site_general_stat(*uncorrelated, D2, total, polarized=True)
[[0.0308642 0.       ]
 [0.        0.0308642]]
>>> _two_site_general_stat(*uncorrelated_biallelic, D2, total, polarized=True)
[[0.0625 0.    ]
 [0.     0.0625]]
>>> _two_site_general_stat(*repulsion_biallelic, D2, total, polarized=True)
[[0.0625 0.0625]
 [0.0625 0.0625]]

>>> _two_site_general_stat(*correlated, D2, total, polarized=False)
[[0.0238446 0.0238446]
 [0.0238446 0.0238446]]
>>> _two_site_general_stat(*correlated_biallelic, D2, total, polarized=False)
[[0.0625 0.0625]
 [0.0625 0.0625]]
>>> _two_site_general_stat(*uncorrelated, D2, total, polarized=False)
[[0.02469136 0.        ]
 [0.         0.02469136]]
>>> _two_site_general_stat(*uncorrelated_biallelic, D2, total, polarized=False)
[[0.0625 0.    ]
 [0.     0.0625]]
>>> _two_site_general_stat(*repulsion_biallelic, D2, total, polarized=False)
[[0.0625 0.0625]
 [0.0625 0.0625]]

>>> _two_site_general_stat(*correlated, Dz, total, polarized=True)
[[ 0.01105014 -0.00556318]
 [-0.00556318  0.00419143]]
>>> _two_site_general_stat(*correlated_biallelic, Dz, total, polarized=True)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*uncorrelated, Dz, total, polarized=True)
[[0.00617284 0.        ]
 [0.         0.00617284]]
>>> _two_site_general_stat(*uncorrelated_biallelic, Dz, total, polarized=True)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*repulsion_biallelic, Dz, total, polarized=True)
[[0. 0.]
 [0. 0.]]

>>> _two_site_general_stat(*correlated, Dz, total, polarized=False)
[[0.00338702 0.00338702]
 [0.00338702 0.00338702]]
>>> _two_site_general_stat(*correlated_biallelic, Dz, total, polarized=False)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*uncorrelated, Dz, total, polarized=False)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*uncorrelated_biallelic, Dz, total, polarized=False)
[[0. 0.]
 [0. 0.]]
>>> _two_site_general_stat(*repulsion_biallelic, Dz, total, polarized=False)
[[0. 0.]
 [0. 0.]]

>>> _two_site_general_stat(*correlated, pi2, total, polarized=True)
[[0.04404816 0.0492303 ]
 [0.0492303  0.0550221 ]]
>>> _two_site_general_stat(*correlated_biallelic, pi2, total, polarized=True)
[[0.0625 0.0625]
 [0.0625 0.0625]]
>>> _two_site_general_stat(*uncorrelated, pi2, total, polarized=True)
[[0.04938272 0.04938272]
 [0.04938272 0.04938272]]
>>> _two_site_general_stat(*uncorrelated_biallelic, pi2, total, polarized=True)
[[0.0625 0.0625]
 [0.0625 0.0625]]
>>> _two_site_general_stat(*repulsion_biallelic, pi2, total, polarized=True)
[[0.0625 0.0625]
 [0.0625 0.0625]]

>>> _two_site_general_stat(*correlated, pi2, total, polarized=False)
[[0.04579248 0.04579248]
 [0.04579248 0.04579248]]
>>> _two_site_general_stat(*correlated_biallelic, pi2, total, polarized=False)
[[0.0625 0.0625]
 [0.0625 0.0625]]
>>> _two_site_general_stat(*uncorrelated, pi2, total, polarized=False)
[[0.04938272 0.04938272]
 [0.04938272 0.04938272]]
>>> _two_site_general_stat(*uncorrelated_biallelic, pi2, total, polarized=False)
[[0.0625 0.0625]
 [0.0625 0.0625]]
>>> _two_site_general_stat(*repulsion_biallelic, pi2, total, polarized=False)
[[0.0625 0.0625]
 [0.0625 0.0625]]

>>> _two_site_general_stat(*correlated, r, hap_weighted, polarized=True)
[[0.66666667 0.44444444]
 [0.44444444 0.77777778]]
>>> _two_site_general_stat(*correlated_biallelic, r, hap_weighted, polarized=True)
[[0.5 0.5]
 [0.5 0.5]]
>>> _two_site_general_stat(*uncorrelated, r, hap_weighted, polarized=True)
[[0.66666667 0.        ]
 [0.         0.66666667]]
>>> _two_site_general_stat(*uncorrelated_biallelic, r, hap_weighted, polarized=True)
[[0.5 0. ]
 [0.  0.5]]
>>> _two_site_general_stat(*repulsion_biallelic, r, hap_weighted, polarized=True)
[[0.5 0. ]
 [0.  0.5]]

>>> _two_site_general_stat(*correlated, r, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]
>>> _two_site_general_stat(*correlated_biallelic, r, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]
>>> _two_site_general_stat(*uncorrelated, r, hap_weighted, polarized=False)
[[1. 0.]
 [0. 1.]]
>>> _two_site_general_stat(*uncorrelated_biallelic, r, hap_weighted, polarized=False)
[[1. 0.]
 [0. 1.]]
>>> _two_site_general_stat(*repulsion_biallelic, r, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]

>>> _two_site_general_stat(*correlated, r2, hap_weighted, polarized=True)
[[0.66666667 0.44444444]
 [0.44444444 0.77777778]]
>>> _two_site_general_stat(*correlated_biallelic, r2, hap_weighted, polarized=True)
[[0.5 0.5]
 [0.5 0.5]]
>>> _two_site_general_stat(*uncorrelated, r2, hap_weighted, polarized=True)
[[0.66666667 0.        ]
 [0.         0.66666667]]
>>> _two_site_general_stat(*uncorrelated_biallelic, r2, hap_weighted, polarized=True)
[[0.5 0. ]
 [0.  0.5]]
>>> _two_site_general_stat(*repulsion_biallelic, r2, hap_weighted, polarized=True)
[[0.5 0. ]
 [0.  0.5]]

>>> _two_site_general_stat(*correlated, r2, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]
>>> _two_site_general_stat(*correlated_biallelic, r2, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]
>>> _two_site_general_stat(*uncorrelated, r2, hap_weighted, polarized=False)
[[1. 0.]
 [0. 1.]]
>>> _two_site_general_stat(*uncorrelated_biallelic, r2, hap_weighted, polarized=False)
[[1. 0.]
 [0. 1.]]
>>> _two_site_general_stat(*repulsion_biallelic, r2, hap_weighted, polarized=False)
[[1. 1.]
 [1. 1.]]

"""
from enum import Enum

import numpy as np


def leave_out_wraparound(n, skip_idx):
    """
    Loop in a circle, skipping the given index.

    >>> list(leave_out_wraparound(5, 0))
    [1, 2, 3, 4]
    >>> list(leave_out_wraparound(5, 1))
    [2, 3, 4, 0]
    >>> list(leave_out_wraparound(5, 2))
    [3, 4, 0, 1]
    >>> list(leave_out_wraparound(0, 0))
    []
    >>> list(leave_out_wraparound(1, 0))
    []
    >>> list(leave_out_wraparound(2, 0))
    [1]
    >>> list(leave_out_wraparound(2, 1))
    [0]

    :param n: length of iterable to iterate over
    :param skip_idx: index to skip when looping

    """
    for i in range(n - 1):
        yield (i + skip_idx + 1) % n


def pairs_with_replacement_idx(n):
    """
    Generate indices for n choose 2 combinations, given n

    >>> from itertools import combinations_with_replacement
    >>> list(pairs_with_replacement_idx(100)) == list(combinations_with_replacement(range(100), 2))
    True
    >>> list(pairs_with_replacement_idx(3))
    [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]

    """
    subloop_start = 0
    for i in range(n):
        for j in range(subloop_start, n):
            yield i, j
        subloop_start += 1


def get_state(ts, state, num_states):
    """
    :param ts: tree sequence object

    """
    for tree in ts.trees():
        for site in tree.sites():
            current_state = 0
            for mutation in site.mutations:
                # TODO: if we haven't encountered this allele before, otherwise current_state = index_of_allele
                current_state += 1
                for sample_id in tree.samples(mutation.node):
                    state[site.id, sample_id] = current_state
            num_states[site.id] = current_state + 1


def get_allele_weights(a_state, b_state, num_samples, haplotype_states):
    """
    >>> haplotype_states = np.zeros((2, 3))
    >>> get_allele_weights([0, 0, 1, 0, 0], [1, 2, 0, 0, 1], 5, haplotype_states)
    >>> haplotype_states
    array([[1., 2., 1.],
           [1., 0., 0.]])

    """
    for i in range(num_samples):
        haplotype_states[a_state[i], b_state[i]] += 1


def two_site_general_stat(ts, summary_func, norm_method, polarized=False, debug=False):
    """
    Wrapper for python prototyping, inline in final version
    """
    state = np.zeros((ts.num_sites, ts.num_samples), dtype=int)
    num_states = np.zeros(ts.num_sites, dtype=int)
    get_state(ts, state, num_states)
    return _two_site_general_stat(state, num_states, summary_func, norm_method, polarized=polarized, debug=debug)


class NormMethod(Enum):
    """
    Method for norming the output statistic

    TOTAL: Divide the statistic by the total number of haplotypes
    HAP_WEIGHTED: Weight each allele's statistic by the proportion of the haplotype present
    """

    TOTAL = "total"
    HAP_WEIGHTED = "hap_weighted"


def _two_site_general_stat(
    state, num_states, summary_func, norm_method, polarized=False, debug=False
):  # pylint:disable=too-many-arguments,too-many-locals,too-many-branches
    norm = NormMethod(norm_method)
    result = np.zeros((len(state), len(state)))
    for a_site_id, b_site_id in pairs_with_replacement_idx(len(state)):
        # TODO: compute_general_two_site_stat_result
        a_dim = num_states[a_site_id]
        b_dim = num_states[b_site_id]

        haplotype_states = np.zeros((a_dim, b_dim), dtype=int)
        get_allele_weights(state[a_site_id], state[b_site_id], state.shape[1], haplotype_states)

        stat_result = 0  # TODO: what is a sensible default here?
        if debug:
            print(f"a: {state[a_site_id]} b: {state[b_site_id]}")
            print(haplotype_states)
            if norm is NormMethod.HAP_WEIGHTED:
                print("wAB\twAb\twaB\tstat\thap_prop")
            else:
                print("wAB\twAb\twaB")
        for a_idx in range(1 if polarized else 0, a_dim):
            for b_idx in range(1 if polarized else 0, b_dim):
                w_AB = haplotype_states[a_idx, b_idx]
                w_Ab = sum(haplotype_states[a_idx, idx] for idx in leave_out_wraparound(b_dim, b_idx))
                w_aB = sum(haplotype_states[idx, b_idx] for idx in leave_out_wraparound(a_dim, a_idx))

                if norm is NormMethod.TOTAL:
                    stat_result += summary_func(w_AB, w_Ab, w_aB, state.shape[1])
                    if debug:
                        print(w_AB, w_Ab, w_aB, summary_func(w_AB, w_Ab, w_aB, state.shape[1]))
                elif norm is NormMethod.HAP_WEIGHTED:
                    haplotype_prop = haplotype_states[a_idx, b_idx] / haplotype_states.sum().sum()

                    stat_result += haplotype_prop * summary_func(w_AB, w_Ab, w_aB, state.shape[1])
                    if debug:
                        print(
                            w_AB, w_Ab, w_aB, summary_func(w_AB, w_Ab, w_aB, state.shape[1]), haplotype_prop, sep="\t"
                        )
                else:
                    raise ValueError(f"unhandled norm method type: {norm}")
        if debug:
            print("--------------------------------------------------------")

        if norm is NormMethod.TOTAL:
            stat_result /= (a_dim - (1 if polarized else 0)) * (b_dim - (1 if polarized else 0))
        result[a_site_id, b_site_id] = stat_result

        # TODO: how to do this in C?
        # TODO: actually we won't, this will be handled by downstream python machinery, so this method is fine.
        # reflect upper triangle to lower triangle
        lower_triangle_idx = np.tril_indices(len(result), k=-1)
        result[lower_triangle_idx] = result.T[lower_triangle_idx]

    return result
