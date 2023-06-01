from enum import Enum

import numpy as np


def pairs_with_replacement(n):
    """
    Generate indices for n choose 2 combinations, given n

    >>> from itertools import combinations_with_replacement
    >>> list(pairs_with_replacement(100)) == list(combinations_with_replacement(range(100), 2))
    True
    >>> list(pairs_with_replacement(3))
    [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]

    """
    subloop_start = 0
    for i in range(n):
        for j in range(subloop_start, n):
            yield i, j
        subloop_start += 1


class NormMethod(Enum):
    """
    Method for norming the output statistic

    TOTAL: Divide the statistic by the total number of haplotypes
    HAP_WEIGHTED: Weight each allele's statistic by the proportion of the haplotype present
    AF_WEIGHTED: Weight each allele's statistic by the product of the allele frequencies
    """

    TOTAL = "total"
    HAP_WEIGHTED = "hap_weighted"
    AF_WEIGHTED = "af_weighted"


def compute_stat_and_weights(hap_mat, summary_func, polarized, norm_method):
    hap_mat = np.asarray(hap_mat)

    # number of samples
    n = hap_mat.sum()
    # number of B alleles, number of A alleles
    _, n_b, n_a = hap_mat.shape

    weights = np.zeros(hap_mat.shape)
    stats = np.zeros(hap_mat.shape)
    for a_idx in range(1 if polarized else 0, n_a):
        for b_idx in range(1 if polarized else 0, n_b):
            # NB: the A, B indices are correct. We are representing A as columns and B as rows
            #     so the indexing looks funny, but it's not
            w_AB = hap_mat[:, b_idx, a_idx]
            w_Ab = hap_mat[:, :, a_idx].sum() - w_AB
            w_aB = hap_mat[:, b_idx, :].sum() - w_AB
            stats[:, b_idx, a_idx] = summary_func(w_AB, w_Ab, w_aB, n)

            # create weights matrix
            if norm_method is NormMethod.HAP_WEIGHTED:
                hap_freq = hap_mat / n
                weights[:, b_idx, a_idx] = hap_freq[:, b_idx, a_idx]
            elif norm_method is NormMethod.AF_WEIGHTED:
                p_a = hap_mat.sum(1) / n
                p_b = hap_mat.sum(2) / n
                weights[:, b_idx, a_idx] = p_a[a_idx] * p_b[b_idx]
            elif norm_method is NormMethod.TOTAL:
                weights[:, b_idx, a_idx] = 1 / ((n_a - (1 if polarized else 0)) * (n_b - (1 if polarized else 0)))

    return stats, weights


def compute_two_site_general_stat(
    state, func, polarized, norm_method, sample_sets, windows, debug=False
):  # pylint:disable=too-many-arguments,too-many-locals
    state = np.asarray(state)
    norm = NormMethod(norm_method)

    n_sample_sets, n_sites = state.shape[:2]
    result = np.zeros((n_sample_sets, n_sites, n_sites))

    for l_idx, r_idx in pairs_with_replacement(n_sites):
        import ipdb; ipdb.set_trace()
        hap_mat = np.zeros((n_sample_sets, len(np.unique(right)), len(np.unique(left))))
        for sample_set_idx, sample_set in enumerate(sample_sets):
            # NB: the left, right / A, B indices are correct. We are representing
            #     A as columns and B as rows, so the indexing looks funny, but it's not
            for a_idx, b_idx in zip(left[sample_set], right[sample_set]):
                hap_mat[sample_set_idx, b_idx, a_idx] += 1
            stats, weights = compute_stat_and_weights(hap_mat, func, polarized, norm)
            if debug:
                print("hap_mat", hap_mat, "stats", stats, "weights", weights, "============", sep="\n")
            result[l_idx, r_idx] = (stats * weights).sum()

    tri_idx = np.tril_indices(len(result), k=-1)
    result[tri_idx] = result.T[tri_idx]
    return result


def get_state(ts, state_dim):
    state = np.zeros((state_dim, ts.num_sites, ts.num_samples), dtype=int)
    for tree in ts.trees():
        for site in tree.sites():
            states = [site.ancestral_state]
            current_state = 0
            for mutation in site.mutations:
                # if we've seen the state, use the index as the enumerated
                # state value, otherwise, add another state value
                if mutation.derived_state in states:
                    current_state = states.index(mutation.derived_state)
                else:
                    states.append(mutation.derived_state)
                    current_state = len(states) - 1

                for sample_id in tree.samples(mutation.node):
                    state[:, site.id, sample_id] = current_state
    return state


def two_site_general_stat(
    ts, summary_func, norm_method, polarized, sample_sets=None, windows=None, debug=False
):  # pylint:disable=too-many-arguments
    sample_sets = sample_sets or np.asarray([ts.samples()])
    assert all(s in ts.samples() for s_set in sample_sets for s in s_set), "bad samples in sample set"
    windows = windows or []  # TODO: implement windows
    state = get_state(ts, len(sample_sets))
    return compute_two_site_general_stat(state, summary_func, polarized, norm_method, sample_sets, windows, debug)
