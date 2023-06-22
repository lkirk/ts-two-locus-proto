project('two-locus-proto', 'c', version: '0.0.1')

c_args = [
    '-g', '-Wall', '-Wextra', '-Werror', '-Wpedantic', '-W',
    '-Wmissing-prototypes',  '-Wstrict-prototypes',
    '-Wconversion', '-Wshadow', '-Wpointer-arith', '-Wcast-align',
    '-Wcast-qual', '-Wwrite-strings', '-Wnested-externs',
    '-fshort-enums', '-fno-common'
]

inc = include_directories(['.'])
deps = [dependency('tskit')]
sources = ['prototype.c']

executable('two-locus-prototype', sources + ['main.c'], c_args: c_args, dependencies: deps, include_directories: inc)

prototype_lib = static_library('ts-prototype', sources: sources, dependencies: deps)

cunit_dep = dependency('cunit')
test_lib = static_library('testlib', sources: ['testlib.c'], dependencies: [cunit_dep, deps])
test_prototype = executable('test_prototype', ['test_prototype.c'],
                            link_with: [prototype_lib, test_lib], c_args: c_args, dependencies: [cunit_dep, deps])
test('prototype', test_prototype)
