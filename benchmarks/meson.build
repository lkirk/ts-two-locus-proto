project('two-locus-benchmark', 'c', version: '0.0.1')

c_args = [
    '-g', '-Wall', '-Wextra', '-Werror', '-Wpedantic', '-W',
    '-Wmissing-prototypes',  '-Wstrict-prototypes',
    '-Wconversion', '-Wshadow', '-Wpointer-arith', '-Wcast-align',
    '-Wcast-qual', '-Wwrite-strings', '-Wnested-externs',
    '-fshort-enums', '-fno-common'
]

tskit_project = subproject('tskit')
tskit_dep = tskit_project.get_variable('tskit_dep')

executable('two-locus-benchmark', ['main.c'], c_args: c_args, dependencies: tskit_dep)
