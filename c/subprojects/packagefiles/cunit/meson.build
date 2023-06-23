# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2020, MayaData Inc. All rights reserved.
# Copyright (c) 2020, DataCore Software Corporation. All rights reserved.

project('cunit', 'C',
	version: '2.1-3',
	default_options: [
		'buildtype=release',
		'default_library=static',
		'b_vscrt=none',
		'warning_level=3'
		],
	meson_version: '>= 0.55.0'
)

cc = meson.get_compiler('c')

includes = include_directories('CUnit', 'CUnit/Headers')

subdir('CUnit')

cunit_dep = declare_dependency(
  include_directories : include_directories('.'),
  link_with : static_lib,
)
