# 2-site stats

Prototype code for 2-site statistics in tskit

## Installation

To use this code/notebook, you'll need to recreate the conda environment like so:

```
conda env create -n <name-of-your-env> --file dev-conda-env.yml
pip install -e .
```

To build the c code:

```
cd c
meson setup build
```

## Running Tests

There are two test suites: a python and a c test suite for running the python
prototype and c prototype tests.

To run the c tests:

```
cd c
ninja -C build test
```

To run the python tests:

```
conda activate <name-of-your-env>
pytest -v test
```

These tests are also run in github actions. The c code coverage can be viewed [here](https://lkirk.github.io/ts-two-locus-proto/)

If you'd like to generate code coverage locally:

```
conda activate <name-of-your-env>
meson setup -Db_coverage=true build  # you might want to remove previous build dir
ninja -C build test coverage
```
