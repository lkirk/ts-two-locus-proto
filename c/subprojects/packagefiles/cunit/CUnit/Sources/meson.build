# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2020, MayaData Inc. All rights reserved.
# Copyright (c) 2020, DataCore Software Corporation. All rights reserved.

cflags = [
    '-D_CRT_SECURE_NO_WARNINGS',
    '-DCU_BUILD_DLL'
]

sources = files(
  'Automated/Automated.c',
  'Basic/Basic.c',
  'Console/Console.c',
  'Framework/CUError.c',
  'Framework/MyMem.c',
  'Framework/TestDB.c',
  'Framework/TestRun.c',
  'Framework/Util.c',
  # 'Win/Win.c',
)

if cc.has_argument('-Wno-unused-but-set-variable')
  cflags += '-Wno-unused-but-set-variable'
endif

if cc.has_argument('-Wno-clobbered')
  cflags += '-Wno-clobbered'
endif

if cc.has_argument('-Wno-format-overflow')
  cflags += '-Wno-format-overflow'
endif

# build static lib
static_lib = static_library('cunit',
		sources,
		c_args: cflags,
		include_directories: [ includes ],
		install: true)

test_cflags = [
    '-DCUNIT_BUILD_TESTS',
    '-DMEMTRACE',
    '-DCUNIT_DO_NOT_DEFINE_UNLESS_BUILDING_TESTS'
]

# build unit tests
test_unit = executable('test_cunit',
    [ 'Test/test_cunit.c', sources ],
		c_args: [ cflags, test_cflags ],
		include_directories: [ 'Test', includes ],
		install: true
)
