def get_paper_ex_ts():
    import io  # pylint:disable=import-outside-toplevel
    import tskit  # pylint:disable=import-outside-toplevel

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

    return tskit.load_text(
        nodes=io.StringIO(nodes),
        edges=io.StringIO(edges),
        sites=io.StringIO(sites),
        individuals=io.StringIO(individuals),
        mutations=io.StringIO(mutations),
        strict=False,
    )
