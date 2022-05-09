gauss2dfit_sources = [
    'centroidparameters.cc',
    'channel.cc',
    'ellipseparameters.cc',
    'ellipticalcomponent.cc',
    'fractionalintegralmodel.cc',
    'gaussiancomponent.cc',
    'integralmodel.cc',
    'linearintegralmodel.cc',
    'psfmodel.cc',
    'source.cc',
    'transforms.cc',
]
gauss2dfitlib = shared_library('gauss2dfit',
                                gauss2dfit_sources,
                                include_directories : inc,
                                dependencies : depends,
                                install : true,
                                soversion : '0.1.0')
