# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2020, MayaData Inc. All rights reserved.
# Copyright (c) 2020, DataCore Software Corporation. All rights reserved.

headers = [
  'Automated.h',
  'Basic.h',
  'Console.h',
  'CUCurses.h',
  'CUError.h',
  'MyMem.h',
  'TestDB.h',
  'TestRun.h',
  'Util.h',
  'wxWidget.h',
]

foreach h : headers
  configure_file(input : 'Headers' / h, output : h, copy : true)
  install_headers('Headers' / h, subdir : 'CUnit')
endforeach

conf_data = configuration_data()
ver = meson.project_version().split('-')
conf_data.set('VERSION', ver.get(0))
conf_data.set('RELEASE', ver.get(1))

cunit_h = configure_file(input : 'Headers/CUnit.h.in',
               output : 'CUnit.h',
               configuration : conf_data)

install_headers(cunit_h, subdir : 'CUnit')

subdir('Sources')