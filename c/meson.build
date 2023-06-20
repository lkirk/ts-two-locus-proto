project('two-locus-proto', 'c', version: '0.0.1')

c_args = [
    '-g', '-Wall', '-Wextra', '-Werror', '-Wpedantic', '-W',
    '-Wmissing-prototypes',  '-Wstrict-prototypes',
    '-Wconversion', '-Wshadow', '-Wpointer-arith', '-Wcast-align',
    '-Wcast-qual', '-Wwrite-strings', '-Wnested-externs',
    '-fshort-enums', '-fno-common'
]

executable('two-locus-debug', ['debug.c'], dependencies: [dependency('tskit')])
executable('two-locus-prototype', ['prototype.c'], dependencies: [dependency('tskit')])
