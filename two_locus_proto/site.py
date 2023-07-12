from enum import Enum
from itertools import combinations_with_replacement

import numpy as np


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


def get_state(ts):
    state = np.zeros((ts.num_sites, ts.num_samples), dtype=int)
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
                    state[site.id, sample_id] = current_state
    return state


def compute_stat_and_weights(hap_mat, summary_func, polarized, norm_method, left_site, right_site, print_weights):
    hap_mat = np.asarray(hap_mat)

    # number of samples
    n = hap_mat.sum()
    # number of B alleles, number of A alleles
    n_a, n_b = hap_mat.shape

    weights = np.zeros(hap_mat.shape)
    stats = np.zeros(hap_mat.shape)
    for a_idx in range(1 if polarized else 0, n_a):
        for b_idx in range(1 if polarized else 0, n_b):
            w_AB = hap_mat[a_idx, b_idx]
            w_Ab = hap_mat[a_idx, :].sum() - w_AB
            w_aB = hap_mat[:, b_idx].sum() - w_AB

            stats[a_idx, b_idx] = summary_func(w_AB, w_Ab, w_aB, n)

            if print_weights:
                print(left_site, right_site, a_idx, b_idx, int(w_AB), int(w_Ab), int(w_aB), int(n), sep="\t")
            # create weights matrix
            if norm_method is NormMethod.HAP_WEIGHTED:
                hap_freq = hap_mat / n
                weights[a_idx, b_idx] = hap_freq[a_idx, b_idx]
            elif norm_method is NormMethod.AF_WEIGHTED:
                p_a = hap_mat.sum(1) / n
                p_b = hap_mat.sum(0) / n
                weights[a_idx, b_idx] = p_a[a_idx] * p_b[b_idx]
            elif norm_method is NormMethod.TOTAL:
                weights[a_idx, b_idx] = 1 / ((n_a - (1 if polarized else 0)) * (n_b - (1 if polarized else 0)))

    return stats, weights


def compute_two_site_general_stat(state, func, polarized, norm_method, debug=False, print_weights=False):
    state = np.asarray(state)
    norm = NormMethod(norm_method)
    # the length of state is the number of sites
    result = np.zeros((len(state), len(state)))
    for (l_idx, left), (r_idx, right) in combinations_with_replacement(enumerate(state), 2):
        hap_mat = np.zeros((np.max(left) + 1, np.max(right) + 1))
        for A_i, B_i in zip(left, right):
            hap_mat[A_i, B_i] += 1
        stats, weights = compute_stat_and_weights(hap_mat, func, polarized, norm, l_idx, r_idx, print_weights)
        if debug:
            print("hap_mat", hap_mat, "stats", stats, "weights", weights, "============", sep="\n")
        result[l_idx, r_idx] = (stats * weights).sum()

    if print_weights:
        np.set_printoptions(precision=15)
        print("result", ",".join(v.astype(str) for v in result[np.triu_indices(len(result))]))
    tri_idx = np.tril_indices(len(result), k=-1)
    result[tri_idx] = result.T[tri_idx]
    return result


def two_site_general_stat(ts, summary_func, norm_method, polarized, debug=False, print_weights=False):
    state = get_state(ts)
    return compute_two_site_general_stat(state, summary_func, polarized, norm_method, debug, print_weights)
