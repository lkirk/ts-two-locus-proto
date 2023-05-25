import io

import pytest
import tskit


@pytest.fixture
def tree_sequence(tree_dict):
    return tskit.load_text(strict=False, **{k: io.StringIO(v) for k, v in tree_dict.items()})
