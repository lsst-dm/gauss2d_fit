test_sources = ['test.cc']
test = executable('gauss2dfit_test',
                  test_sources,
                  dependencies: depends,
                  include_directories : inc,
                  link_with : gauss2dfitlib,
                  )
test('Tests', test)