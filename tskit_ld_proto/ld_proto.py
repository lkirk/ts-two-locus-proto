import io
from enum import Enum

import numpy as np
import tskit

# Data taken from the tests: https://github.com/tskit-dev/tskit/blob/61a844a/c/tests/testlib.c#L55-L96

nodes = """\
is_sample time population individual
1  0       -1   0
1  0       -1   0
1  0       -1   1
1  0       -1   1
0  0.071   -1   -1
0  0.090   -1   -1
0  0.170   -1   -1
0  0.202   -1   -1
0  0.253   -1   -1
"""

edges = """\
left   right   parent  child
2 10 4 2
2 10 4 3
0 10 5 1
0 2  5 3
2 10 5 4
0 7  6 0,5
7 10 7 0,5
0 2  8 2,6
"""

sites = """\
position ancestral_state
1      0
4.5    0
8.5    0
"""

mutations = """\
site node derived_state
0      2   1
1      0   1
2      5   1
"""

individuals = """\
flags  location   parents
0      0.2,1.5    -1,-1
0      0.0,0.0    -1,-1
"""

paper_ex_ts = tskit.load_text(
    nodes=io.StringIO(nodes),
    edges=io.StringIO(edges),
    sites=io.StringIO(sites),
    individuals=io.StringIO(individuals),
    mutations=io.StringIO(mutations),
    strict=False,
)


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


def compute_D(w_AB, w_Ab, w_aB, n):
    """
    >>> compute_D(3, 0, 0, 4)
    0.1875

    >>> compute_D(2, 1, 1, 4)
    -0.0625

    """
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    return p_AB - (p_A * p_B)


def compute_r2(w_AB, w_Ab, w_aB, n):
    """
    >>> compute_r2(3, 0, 0, 4)
    1.0

    >>> compute_r2(2, 1, 1, 4)
    0.1111111111111111

    """
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    D = p_AB - (p_A * p_B)
    denom = p_A * p_B * (1 - p_A) * (1 - p_B)

    if denom == 0 and D == 0:
        return np.nan

    return (D * D) / denom


def two_site_general_stat(ts, summary_func, norm_method):
    """
    Wrapper for python prototyping, inline in final version
    """
    state = np.zeros((ts.num_sites, ts.num_samples), dtype=int)
    num_states = np.zeros(ts.num_sites, dtype=int)
    get_state(ts, state, num_states)
    return _two_site_general_stat(state, num_states, summary_func, norm_method)


class NormMethod(Enum):
    """
    Method for norming the output statistic

    TOTAL: Divide the statistic by the total number of haplotypes
    HAP_WEIGHTED: Weight each allele's statistic by the proportion of the haplotype present
    """

    TOTAL = "total"
    HAP_WEIGHTED = "hap_weighted"


def _two_site_general_stat(state, num_states, summary_func, norm_method, polarized=False, debug=False):
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
