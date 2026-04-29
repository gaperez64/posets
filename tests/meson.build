tests_exe = executable ('tests', 'tests.cc',
                        dependencies : [posets_dep ],
                        c_args : '-O0' )  # Deactivate optimization, otherwise it takes ages.

kdtests_exe = executable ('kdtests', 'kdtests.cc',
                          dependencies : [posets_dep ])

sttests_exe = executable ('sttests', 'sttests.cc',
                          dependencies : [posets_dep ])

csttests_exe = executable ('csttests', 'csttests.cc',
                           dependencies : [posets_dep ])

simdbm_exe = executable ('simd-bm', 'simd-bm.cc',
                         dependencies : [posets_dep ])

downsetbm_exe = executable ('downset-bm', 'downset-bm.cc',
                            dependencies : [posets_dep ])

lcrstests_exe = executable ('lcrstests', 'lcrstests.cc',
                            dependencies : [posets_dep ])

test('antichains/vectors random', downsetbm_exe,
     args : ['all', 'all', '--params=build=5,query=5,transfer=5,intersection=5,union=5' ])
test('antichains/vectors implementations', tests_exe, args : ['all', 'all'])
test('antichains/vectors kdtests', kdtests_exe)
test('antichains/vectors sttests', sttests_exe)
test('antichains/vectors csttests', csttests_exe)
test('antichains/vectors lcrstests', lcrstests_exe)
