import io

import msprime
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
    strict=False
)

def leave_out_wraparound(num_states, skip_state):
    for i in range(num_states - 1):
        yield (i + skip_state + 1) % num_states


def pairs_with_replacement_idx(num):
    l = num
    s = 0
    for i in range(l):
        for j in range(s, l):
            yield i, j
        s += 1


from itertools import combinations_with_replacement
assert all([l == r for l, r in zip(combinations_with_replacement(range(100), 2), pairs_with_replacement_idx(100))]), 'test fails'


def get_state(ts):
    state = np.zeros((ts.num_sites, ts.num_samples), dtype=int)
    num_states = np.zeros(ts.num_sites, dtype=int)
    for tree in ts.trees():
        for site in tree.sites():
            current_state = 0
            for mutation in site.mutations:
                # TODO: if we haven't encountered this allele before, otherwise current_state = index_of_allele
                current_state += 1
                for sample_id in tree.samples(mutation.node):
                    state[site.id, sample_id] = current_state
            num_states[site.id] = current_state
    return num_states, state


def get_allele_weights(state, num_states):
    for (A_site_idx, A_samples), (B_site_idx, B_samples) in list(combinations_with_replacement(enumerate(state), 2)):
        # print(f'{A_site_idx},{B_site_idx},{A_samples},{B_samples}')
        A_dim = num_states[A_site_idx]
        B_dim = num_states[B_site_idx]

        haplotype_counts = np.zeros((A_dim + 1, B_dim + 1), dtype=int)
        for A_state, B_state in zip(A_samples, B_samples):
            haplotype_counts[A_state, B_state] += 1

        # print(haplotype_counts)
        for A in range(A_dim):
            for B in range(B_dim):
                w_AB = haplotype_counts[A, B]
                w_Ab = sum([haplotype_counts[A, idx] for idx in leave_out_wraparound(B_dim + 1, B)])
                w_aB = sum([haplotype_counts[idx, B] for idx in leave_out_wraparound(A_dim + 1, A)])
                yield A_site_idx, B_site_idx, (w_AB, w_Ab, w_aB)


# TODO: why wrong?
# def compute_D(w_AB, w_Ab, w_aB, n):
#     p_AB = w_AB / float(n)
#     p_Ab = w_Ab / float(n)
#     p_aB = w_aB / float(n)
#     p_ab = 1 - p_AB
#     return p_AB * p_ab - p_Ab * p_aB


# TODO: why wrong?
# def compute_D(w_AB, w_Ab, w_aB, n):
#     p_AB = w_AB / float(n)
#     p_Ab = w_Ab / float(n)
#     p_aB = w_aB / float(n)
#     p_ab = 1 - p_AB
#     return p_AB * p_ab - p_Ab * p_aB


def compute_D(w_AB, w_Ab, w_aB, n):
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    return p_AB - (p_A * p_B)


def compute_r2(w_AB, w_Ab, w_aB, n):
    p_AB = w_AB / float(n)
    p_Ab = w_Ab / float(n)
    p_aB = w_aB / float(n)
    p_ab = 1 - p_AB

    p_A = p_AB + p_Ab
    p_B = p_AB + p_aB

    D = p_AB - (p_A * p_B)

    # D = (p_AB * p_ab) - (p_Ab * p_aB)

    denom = p_A * p_B * (1 - p_A) * (1 - p_B)

    if denom == 0 and D == 0:
        return np.nan
    else:
        return (D * D) / denom

def compute_stat_matrix(ts, summary_func):
    num_states, state = get_state(ts)
    matrix = np.zeros((len(state), len(state)))
    for i, j, w in get_allele_weights(state, num_states):
        matrix[i, j] = summary_func(*w, state.shape[1]) 
    lower_triangle_idx = np.tril_indices(len(matrix), k=-1)
    matrix[lower_triangle_idx] = matrix.T[lower_triangle_idx]
    return matrix
