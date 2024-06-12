import lsst.gauss2d as g2d
import lsst.gauss2d.fit as g2f
import numpy as np
import pytest


@pytest.fixture(scope='module')
def channels():
    return tuple(g2f.Channel(x) for x in 'rgb')


@pytest.fixture(scope='module')
def data(channels):
    n_x, n_y = 5, 5
    return g2f.DataD([
        g2f.ObservationD(
            g2d.ImageD(np.zeros((n_y, n_x))),
            g2d.ImageD(np.full((n_y, n_x), 1e3)),
            g2d.ImageB(np.ones((n_y, n_x))),
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
                    [(g2f.Channel.NONE, g2f.IntegralParameterD(1.))]
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
def interptype_sersic_default():
    return g2f.SersicMixComponentIndexParameterD().interptype


@pytest.fixture(scope='module')
def sources(channels, interptype_sersic_default):
    last = None
    n_sources, n_components = 2, 2
    sources = [None] * n_sources
    n_comp_max = n_components - 1
    log10 = g2f.Log10TransformD()
    for i in range(n_sources):
        components = [None] * n_components
        for c in range(n_components):
            is_last = c == n_comp_max
            last = g2f.FractionalIntegralModel(
                [
                    (channel, g2f.ProperFractionParameterD(
                        (is_last == 1) or 0.25 * (1 + idx_channel),
                        fixed=is_last
                    ))
                    for idx_channel, channel in enumerate(channels)
                ],
                g2f.LinearIntegralModel([
                    (channel, g2f.IntegralParameterD(0.5 + 0.5 * (i + 1)))
                    for channel in channels
                ]) if (c == 0) else last,
                is_last,
            )
            components[c] = g2f.GaussianComponent(
                centroid=None,
                ellipse=g2f.GaussianParametricEllipse(1., 1., 0.),
                integral=last,
            )
        components.append(
            g2f.SersicMixComponent(
                centroid=None,
                ellipse=g2f.SersicParametricEllipse(
                    size_x=g2f.ReffXParameterD(2.6, transform=log10, fixed=False),
                    size_y=g2f.ReffYParameterD(3.5, transform=log10, fixed=False),
                    rho=g2f.RhoParameterD(
                        -0.1,
                        fixed=False,
                        transform=g2f.LogitLimitedTransformD(
                            limits=g2f.LimitsD(min=-0.999, max=0.999, name="ref_logit_sersic[0.5, 6.0]"),
                        ),
                    )
                ),
                integral=g2f.LinearIntegralModel([
                    (
                        channel,
                        g2f.IntegralParameterD(1.0, transform=log10, label=f"Sersic {channel}-band")
                    )
                    for channel in channels
                ]),
                # Linear interpolation fails at exact knot values
                # Adding a small offset solves the problem
                sersicindex=g2f.SersicMixComponentIndexParameterD(
                    value=1.0 + 5e-4 * (interptype_sersic_default == g2f.InterpType.linear),
                    fixed=False,
                    transform=g2f.LogitLimitedTransformD(
                        limits=g2f.LimitsD(min=0.5, max=6.0, name="ref_logit_sersic[0.5, 6.0]"),
                    ),
                )
            ),
        )
        sources[i] = g2f.Source(components)

    return sources


@pytest.fixture(scope='module')
def model(data, psfmodels, sources):
    model = g2f.ModelD(data, list(psfmodels), list(sources))
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
    assert len(gaussians) == 12
    params = model.parameters()
    assert len(params) == 80
    assert len(set(params)) == 68


def test_model_eval_jacobian(model):
    model.setup_evaluators(g2f.EvaluatorMode.jacobian)
    result = np.array(model.evaluate())
    assert (result == 0).all()
    # TODO: Investigate why this rtol needs to be so high
    errors = model.verify_jacobian(atol=1e-3, rtol=5e-3)
    # All of the IntegralParameters got double-counted - see DM-40674
    # ... but also, linear interpolators will fail the Sersic index params
    assert len(errors) == 6


def test_model_eval_loglike_grad(model):
    model.setup_evaluators(g2f.EvaluatorMode.loglike_grad, print=True)
    result = np.array(model.evaluate())
    assert (result == 0).all()
    dloglike = np.array(model.compute_loglike_grad(True))
    assert (dloglike == 0).all()


def test_model_hessian(model):
    options = g2f.HessianOptions(return_negative=False)
    for idx in range(2):
        hessians = {is_transformed: model.compute_hessian(transformed=is_transformed, options=options)
                    for is_transformed in (False, True)}
        # Note that since there is no noise in the image and none of the param values have changed,
        # the Hessian terms are not all negative - the sign corrections are based on the sign of loglike_grad,
        # and those are all exactly zero
        assert all(np.isfinite(hessian.data).all() for hessian in hessians.values())
        if idx == 0:
            for param in model.parameters():
                # TODO: Remove this when DM-40674 is fixed
                if isinstance(param, g2f.IntegralParameterD):
                    param.fixed = True

    for is_transformed in (False, True):
        hessian_from_jac = model.compute_hessian(transformed=is_transformed, options=None, print=True).data
        # TODO: investigate these unreasonably large differences after DM-40674
        assert np.all(np.isclose(hessians[is_transformed].data, hessian_from_jac, rtol=1e-3, atol=1))
