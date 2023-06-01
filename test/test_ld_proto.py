import numpy as np
import pytest
import tskit

from two_locus_proto.summary_functions import D, D2, r2, r, D_prime, Dz, pi2
from two_locus_proto.site import compute_two_site_general_stat, get_state, two_site_general_stat

TEST_DATA = dict(
    CORRELATED=[[0, 1, 1, 0, 2, 2, 1, 0, 1], [1, 2, 2, 1, 0, 0, 2, 1, 2]],
    CORRELATED_BIALLELIC=[[0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1]],
    UNCORRELATED=[[0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]],
    UNCORRELATED_BIALLELIC=[[0, 0, 0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, 0, 1, 1]],
    REPULSION_BIALLELIC=[[0, 0, 0, 0, 1, 1, 1, 1], [1, 1, 1, 1, 0, 0, 0, 0]],
)

NORM_METHOD = {
    D: "total",
    D_prime: "hap_weighted",
    D2: "total",
    Dz: "total",
    pi2: "total",
    r: "total",
    r2: "hap_weighted",
}

POLARIZATION = {
    D: True,
    D_prime: True,
    D2: False,
    Dz: False,
    pi2: False,
    r: True,
    r2: False,
}

# fmt: off
TEST_CASES = [
    (D, "CORRELATED",
     np.array([[ 0.05555556, -0.01851852],
               [-0.01851852,  0.04320988]])),
    (D, "CORRELATED_BIALLELIC",
     np.array([[0.25, 0.25],
               [0.25, 0.25]])),
    (D, "UNCORRELATED",
     np.array([[0.05555556, 0.        ],
               [0.        , 0.05555556]])),
    (D, "UNCORRELATED_BIALLELIC",
     np.array([[0.25, 0.  ],
               [0.  , 0.25]])),
    (D, "REPULSION_BIALLELIC",
     np.array([[ 0.25, -0.25],
               [-0.25,  0.25]])),

    (D_prime, "CORRELATED",
     np.array([[0.66666667, 0.44444444],
               [0.44444444, 0.77777778]])),
    (D_prime, "CORRELATED_BIALLELIC",
     np.array([[0.5, 0.5],
               [0.5, 0.5]])),
    (D_prime, "UNCORRELATED",
     np.array([[0.66666667, 0.        ],
               [0.        , 0.66666667]])),
    (D_prime, "UNCORRELATED_BIALLELIC",
     np.array([[0.5, 0. ],
               [0. , 0.5]])),
    (D_prime, "REPULSION_BIALLELIC",
     np.array([[0.5, 0. ],
               [0. , 0.5]])),

    (D2, "CORRELATED",
     np.array([[0.02384460363426, 0.02384460363426],
               [0.02384460363426, 0.02384460363426]])),
    (D2, "CORRELATED_BIALLELIC",
     np.array([[0.0625, 0.0625],
               [0.0625, 0.0625]])),
    (D2, "UNCORRELATED",
     np.array([[0.02469136, 0.        ],
               [0.        , 0.02469136]])),
    (D2, "UNCORRELATED_BIALLELIC",
     np.array([[0.0625, 0.    ],
               [0.    , 0.0625]])),
    (D2, "REPULSION_BIALLELIC",
     np.array([[0.0625, 0.0625],
               [0.0625, 0.0625]])),

    (Dz, "CORRELATED",
     np.array([[0.00338701756168, 0.00338701756168],
               [0.00338701756168, 0.00338701756168]])),
    (Dz, "CORRELATED_BIALLELIC",
     np.array([[0., 0.],
               [0., 0.]])),
    (Dz, "UNCORRELATED",
     np.array([[0., 0.],
               [0., 0.]])),
    (Dz, "UNCORRELATED_BIALLELIC",
     np.array([[0., 0.],
               [0., 0.]])),
    (Dz, "REPULSION_BIALLELIC",
     np.array([[0., 0.],
               [0., 0.]])),

    (pi2, "CORRELATED",
     np.array([[0.04579248, 0.04579248],
               [0.04579248, 0.04579248]])),
    (pi2, "CORRELATED_BIALLELIC",
     np.array([[0.0625, 0.0625],
               [0.0625, 0.0625]])),
    (pi2, "UNCORRELATED",
     np.array([[0.04938272, 0.04938272],
               [0.04938272, 0.04938272]])),
    (pi2, "UNCORRELATED_BIALLELIC",
     np.array([[0.0625, 0.0625],
               [0.0625, 0.0625]])),
    (pi2, "REPULSION_BIALLELIC",
     np.array([[0.0625, 0.0625],
               [0.0625, 0.0625]])),

    (r, "CORRELATED",
     np.array([[ 0.26095428, -0.12212786],
               [-0.12212786,  0.18377223]])),
    (r, "CORRELATED_BIALLELIC",
     np.array([[1., 1.],
               [1., 1.]])),
    (r, "UNCORRELATED",
     np.array([[0.25, 0.  ],
               [0.  , 0.25]])),
    (r, "UNCORRELATED_BIALLELIC",
     np.array([[1., 0.],
               [0., 1.]])),
    (r, "REPULSION_BIALLELIC",
     np.array([[ 1., -1.],
               [-1.,  1.]])),

    (r2, "CORRELATED",
     np.array([[1., 1.],
               [1., 1.]])),
    (r2, "CORRELATED_BIALLELIC",
     np.array([[1., 1.],
               [1., 1.]])),
    (r2, "UNCORRELATED",
     np.array([[1., 0.],
               [0., 1.]])),
    (r2, "UNCORRELATED_BIALLELIC",
     np.array([[1., 0.],
               [0., 1.]])),
    (r2, "REPULSION_BIALLELIC",
     np.array([[1., 1.],
               [1., 1.]])),
]
# fmt: on


@pytest.mark.parametrize(
    "summary_func, test_data, truth", TEST_CASES, ids=[f"{func.__name__}_{data}" for (func, data, _) in TEST_CASES]
)
def test_statistics(summary_func, test_data, truth):
    """
    test the raw statistics functions and methods for computing weights
    do not compute state from tree sequences
    """
    result = compute_two_site_general_stat(
        TEST_DATA[test_data], summary_func, POLARIZATION[summary_func], NORM_METHOD[summary_func]
    )
    np.testing.assert_allclose(result, truth)


TREE_PAPER_EX = dict(
    nodes="""\
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
    """,
    edges="""\
    left   right   parent  child
    2 10 4 2
    2 10 4 3
    0 10 5 1
    0 2  5 3
    2 10 5 4
    0 7  6 0,5
    7 10 7 0,5
    0 2  8 2,6
    """,
    sites="""\
    position ancestral_state
    1      0
    4.5    0
    8.5    0
    """,
    mutations="""\
    site node derived_state
    0      2   1
    1      0   1
    2      5   1
    """,
    individuals="""\
    flags  location   parents
    0      0.2,1.5    -1,-1
    0      0.0,0.0    -1,-1
    """,
)


TREE_BASIC = dict(
    nodes="""\
    is_sample time
    1         0 
    1         0
    1         0
    1         0
    1         0
    0         1
    0         2
    0         3
    0         4
    0         5
    0         6
    0         7
    0         8
    """,
    edges="""\
    left right parent child
    0    100   10      0,1
    0    100   5      2,3
    0    100   8      4,5
    0    100   11     8,10
    100  200   12     0,9
    100  200   9      1,7
    100  200   7      2,6
    100  200   6      3,4
    """,
    sites="""\
    position ancestral_state
    10       A
    30       C
    70       G
    150      G
    170      C
    """,
    mutations="""\
    site node parent derived_state
    0    0    -1     T
    1    8    -1     T
    2    3    -1     T
    3    1    -1     T
    4    4    -1     T
    """,
)

TREE_BACK_MUTATION = {
    **TREE_BASIC,
    **dict(
        sites="""\
        position ancestral_state
        10       A
        30       C
        70       G
        110      A
        150      G
        170      C
        """,
        mutations="""\
        site node parent derived_state
        0    0    -1     T
        1    8    -1     T
        2    3    -1     T
        4    1    -1     T
        5    4    -1     T
        3    9    -1     T
        3    7    5      A
    """,
    ),
}

TREE_BACK_MUTATION_MULTIALLELIC = {
    **TREE_BASIC,
    **dict(
        sites="""\
        position ancestral_state
        10       A
        30       C
        70       G
        110      A
        150      G
        170      C
        """,
        mutations="""\
        site node parent derived_state
        0    0    -1     T
        1    8    -1     T
        2    3    -1     T
        4    1    -1     T
        5    4    -1     T
        1    5    1      A
        3    9    -1     T
        3    7    6      A
        3    2    7      G
    """,
    ),
}


# fmt:off
GET_STATE_CASES = {
    "paper_ex": (TREE_PAPER_EX,
                 np.array([[0, 0, 1, 0],
                           [1, 0, 0, 0],
                           [0, 1, 1, 1]])),
    "tree_basic": (TREE_BASIC,
                   np.array([[1, 0, 0, 0, 0],
                             [0, 0, 1, 1, 1],
                             [0, 0, 0, 1, 0],
                             [0, 1, 0, 0, 0],
                             [0, 0, 0, 0, 1]])),
    "tree_back_mutation": (TREE_BACK_MUTATION,
                           np.array([[1, 0, 0, 0, 0],
                                     [0, 0, 1, 1, 1],
                                     [0, 0, 0, 1, 0],
                                     [0, 1, 0, 0, 0],
                                     [0, 1, 0, 0, 0],
                                     [0, 0, 0, 0, 1]])),
    "tree_back_mutation_multiallelic": (TREE_BACK_MUTATION_MULTIALLELIC,
                           np.array([[1, 0, 0, 0, 0],
                                     [0, 0, 2, 2, 1],
                                     [0, 0, 0, 1, 0],
                                     [0, 1, 2, 0, 0],
                                     [0, 1, 0, 0, 0],
                                     [0, 0, 0, 0, 1]]))
}
# fmt:on


@pytest.mark.parametrize("tree_dict,truth", GET_STATE_CASES.values(), ids=GET_STATE_CASES.keys())
def test_get_state(tree_sequence, truth):
    np.testing.assert_equal(get_state(tree_sequence), truth)


@pytest.mark.parametrize("tree_dict", [TREE_BASIC])
def test_against_tskit(tree_sequence):
    ld_calc = tskit.LdCalculator(tree_sequence).get_r2_matrix()
    result = two_site_general_stat(tree_sequence, r2, NORM_METHOD[r2], POLARIZATION[r2])
    np.testing.assert_allclose(ld_calc, result)
