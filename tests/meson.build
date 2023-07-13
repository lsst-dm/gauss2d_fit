test_names = [
    'centroidparameters',
    'channel',
    'data',
    'fractionalintegralmodel',
    'gaussiancomponent',
    'gaussianparametricellipse',
    'gaussianprior',
    'linearintegralmodel',
    'linearsersicmixinterpolator',
    'model',
    'observation',
    'parametricgaussian1d',
    'psfmodel',
    'sersicmix',
    'sersicmixcomponent',
    'shapeprior',
    'source',
]
foreach test_name : test_names
    test = executable(
        'gauss2dfit_test_' + test_name,
        'test_' + test_name + '.cc',
        dependencies: depends,
        include_directories : inc,
        link_with : gauss2dfitlib,
    )
    test(test_name, test)
endforeach