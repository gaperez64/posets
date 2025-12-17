project ('Posets Library', ['cpp'],
         default_options : ['cpp_std=c++20', 'b_ndebug=if-release'],
         license : 'GPL-3.0')

if not meson.is_subproject()
  extra_args = ['-ftrack-macro-expansion=0',
                '-fno-diagnostics-show-caret',
                '-fconcepts-diagnostics-depth=2']
  add_global_arguments (extra_args, language: 'cpp')
endif

boost_inc = include_directories('lib')

subdir ('include')
subdir ('tests')
