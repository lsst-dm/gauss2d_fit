import gauss2d as g2
import gauss2d.fit as g2f
import numpy as np
import pytest


@pytest.fixture(scope='module')
def channels():
    return tuple(g2f.Channel(x) for x in 'rgb')


@pytest.fixture(scope='module')
def data(channels):
    n_x, n_y = 5, 5
    return g2f.Data([
        g2f.Observation(
            g2.ImageD(np.zeros((n_y, n_x))),
            g2.ImageD(np.full((n_y, n_x), 1e6)),
            g2.ImageB(np.ones((n_y, n_x))),
            channel,
        )
        for channel in channels
    ])


@pytest.fixture(scope='module')
def psfmodels(data):
    psfmodels = tuple(
        g2f.PsfModel([
            g2f.GaussianComponent(
                g2f.GaussianParametricEllipse(1., 1., 0.),
                None,
                g2f.LinearIntegralModel(
                    {g2f.Channel.NONE: g2f.IntegralParameterD(1.)}
                ),
            )
        ])
        for _ in range(len(data))
    )
    for psfmodel in psfmodels:
        for param in psfmodel.parameters():
            param.fixed = True
    return psfmodels


@pytest.fixture(scope='module')
def sources(channels):
    last = None
    n_sources, n_components = 2, 2
    sources = [None]*n_sources
    n_comp_max = n_components - 1
    for i in range(n_sources):
        components = [None]*n_components
        for c in range(n_components):
            is_last = c == n_comp_max
            last = g2f.FractionalIntegralModel(
                {
                    channel: g2f.ProperFractionParameterD(
                        (is_last == 1) or (0.5 + 0.5*(c > 0)),
                        fixed=is_last
                    )
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
    model = g2f.Model(data, list(psfmodels), list(sources))
    with pytest.raises(RuntimeError):
        model.evaluate()
    model.setup_evaluators()
    model.evaluate()
    result = model.evaluate()
    for value in result:
        assert value == 0
    for data, output in zip(model.data, model.outputs):
        data.image.data.flat = output.data.flat
    return model


def test_data(channels, data):
    assert len(data) == len(channels)


def test_model(channels, model):
    assert str(model) != ""
    gaussians = model.gaussians(channels[0])
    assert len(gaussians) == 4
    params = model.parameters()
    assert len(params) == 62
    assert len(set(params)) == 50


def test_model_eval_jacobian(model):
    model.setup_evaluators(g2f.Model.EvaluatorMode.jacobian)
    result = np.array(model.evaluate())
    assert (result == 0).all()


def test_model_eval_loglike_grad(model):
    model.setup_evaluators(g2f.Model.EvaluatorMode.loglike_grad, print=True)
    result = np.array(model.evaluate())
    assert (result == 0).all()
    dloglike = np.array(model.compute_loglike_grad(True))
    assert (dloglike == 0).all()


def test_model_hessian(model):
    hessians = {is_transformed: model.compute_hessian(is_transformed)
                for is_transformed in (False, True)}
    assert all(np.isfinite(hessian.data).all() for hessian in hessians.values())
    for is_transformed in (False, True):
        assert np.all(np.isclose(hessians[is_transformed].data, model.compute_hessian(is_transformed).data))
