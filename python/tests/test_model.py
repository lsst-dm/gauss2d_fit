import gauss2d as g2
import gauss2d.fit as g2f
import pytest


@pytest.fixture(scope='module')
def channels():
    return tuple(g2f.Channel(x) for x in 'rgb')


@pytest.fixture(scope='module')
def data(channels):
    n_x, n_y = 5, 5
    return g2f.Data([
        g2f.Observation(
            g2.ImageD(n_y, n_x),
            g2.ImageD(n_y, n_x),
            g2.ImageB(n_y, n_x),
            channel,
        )
        for channel in channels
    ])


@pytest.fixture(scope='module')
def psfmodels(data):
    return tuple(
        g2f.PsfModel([
            g2f.GaussianComponent(
                g2f.GaussianParametricEllipse(1., 1., 0.),
                None,
                g2f.LinearIntegralModel(
                    {g2f.Channel.NONE: g2f.IntegralParameterD(1.)}
                ),
            )
        ])
        for i in range(len(data))
    )


@pytest.fixture(scope='module')
def sources(channels):
    last = None
    n_sources, n_components = 2, 2
    sources = [None]*n_sources
    n_src_max = n_sources - 1
    for i in range(n_sources):
        components = [None]*n_components
        for c in range(n_components):
            is_last = i == n_src_max
            last = g2f.FractionalIntegralModel(
                {
                    channel: g2f.ProperFractionParameterD((is_last == 1) or (0.5 + 0.5*(c > 0)), fixed=is_last)
                    for channel in channels
                },
                g2f.LinearIntegralModel({
                    channel: g2f.IntegralParameterD(0.5 + 0.5*(i + 1))
                    for channel in channels
                }) if (c == 0) else last,
                is_last,
            )
            components[c] = g2f.GaussianComponent(
                g2f.GaussianParametricEllipse(1., 1., 0.),
                None,
                last,
            )
        sources[i] = g2f.Source(components)

    return sources


@pytest.fixture(scope='module')
def model(data, psfmodels, sources):
    return g2f.Model(data, list(psfmodels), list(sources))


def test_data(channels, data):
    assert len(data) == len(channels)


def test_model(channels, model):
    assert str(model) != ""
    gaussians = model.gaussians(channels[0])
    assert len(gaussians) == 4
    params = model.parameters()
    assert len(params) == 65
    assert len(set(params)) == 54


def test_model_evaluation(model):
    with pytest.raises(RuntimeError):
        model.evaluate()
    model.setup_evaluators()
    result = model.evaluate()
    for value in result:
        assert value == 0
