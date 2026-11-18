import numpy as np
from scipy.integrate import quad

from lsst.gauss2d.fit import IntegralSize
from lsst.gauss2d.fit.utils import (
    gaussian_integrated_flux, gaussian_mixture_integrated_flux, sersic_integrated_flux, sersic_surfbright,
)
from lsst.gauss2d import M_HWHM_SIGMA, M_SIGMA_HWHM

import pytest


@pytest.fixture(scope="module")
def integralsizes() -> tuple[IntegralSize]:
    return (
        IntegralSize(integral=3.246685e-01, sigma=1.444509e+00),
        IntegralSize(integral=4.780321e-01, sigma=8.307142e-01),
        IntegralSize(integral=1.725859e-01, sigma=4.372915e-01),
        IntegralSize(integral=2.471355e-02, sigma=1.865699e-01),
    )


@pytest.fixture(scope="module")
def radii() -> np.ndarray[float]:
    """Return radii for which 0.5<=n<=6.0 profiles vary measurably."""
    return 10 ** np.linspace(-5.1, 0.8, 119)


def test_gaussian_integrated_flux(radii):
    for bad_sigma in (0., -5, -np.inf, np.nan):
        with pytest.raises(ValueError):
            gaussian_integrated_flux(1.0, sigma=bad_sigma)

    # This should be exactly equal, but allow for small roundoff errors
    np.testing.assert_almost_equal(gaussian_integrated_flux(M_HWHM_SIGMA), 0.5, decimal=12)
    assert gaussian_integrated_flux(0) == 0
    assert gaussian_integrated_flux(np.inf) == 1.0

    fluxes = gaussian_integrated_flux(radii)

    # Test monotonic increase in integrated flux
    np.testing.assert_array_less(fluxes[:-1], fluxes[1:])
    np.testing.assert_array_less(fluxes[:-1], 1.0)

    for ratio in (0.90, 1.08):
        fluxes_rescaled = gaussian_integrated_flux(radii*ratio, sigma=ratio)
        np.testing.assert_array_almost_equal(fluxes, fluxes_rescaled, decimal=8)


def test_gaussian_mixture_integrated_flux(integralsizes, radii):
    fluxes = gaussian_mixture_integrated_flux(radii, integralsizes)

    # Test monotonic increase in integrated flux
    np.testing.assert_array_less(fluxes[:-1], fluxes[1:])
    np.testing.assert_array_less(fluxes[:-1], 1.0)


def test_sersic_integrated_flux(radii):
    for bad_sersic_index in (0., -5, -np.inf, np.nan):
        with pytest.raises(ValueError):
            sersic_integrated_flux(1.0, sersic_index=bad_sersic_index)
        for bad_reff in (0., -5, -np.inf, np.nan):
            with pytest.raises(ValueError):
                sersic_integrated_flux(1.0, sersic_index=bad_sersic_index, r_eff=bad_reff)

    # This should be exactly equal, but allow for small roundoff errors
    np.testing.assert_almost_equal(sersic_integrated_flux(1.0, sersic_index=0.5), 0.5, decimal=12)

    for sersic_index in (1e-3, 0.1, 0.5, 1.0, 5.0, 20.):
        assert sersic_integrated_flux(0, sersic_index=sersic_index) == 0
        assert sersic_integrated_flux(np.inf, sersic_index=sersic_index) == 1.0

    for sersic_index in (0.5, 1.0, 2.0, 4.0, 6.0):
        msg = f"{sersic_index=}"
        fluxes = sersic_integrated_flux(radii, sersic_index=sersic_index)

        if sersic_index == 0.5:
            np.testing.assert_array_almost_equal(
                fluxes, gaussian_integrated_flux(radii, sigma=M_SIGMA_HWHM), decimal=8,
            )

        # Test monotonic increase in integrated flux
        np.testing.assert_array_less(fluxes[:-1], fluxes[1:], err_msg=msg)
        # This will work as long as the maximum radius isn't too large
        np.testing.assert_array_less(fluxes, 1.0, err_msg=msg)

        for ratio in (0.92, 1.06):
            fluxes_rescaled = sersic_integrated_flux(radii*ratio, sersic_index=sersic_index, r_eff=ratio)
            np.testing.assert_array_almost_equal(fluxes, fluxes_rescaled, decimal=8, err_msg=msg)


def test_sersic_surfbright(radii):
    for bad_sersic_index in (0., -5, -np.inf, np.nan):
        with pytest.raises(ValueError):
            sersic_surfbright(1.0, sersic_index=bad_sersic_index)
        for bad_reff in (0., -5, -np.inf, np.nan):
            with pytest.raises(ValueError):
                sersic_surfbright(1.0, sersic_index=bad_sersic_index, r_eff=bad_reff)

    radii_start = np.concatenate(([0,], radii[:-1]))

    for sersic_index in (0.5, 1.0, 2.0, 4.0, 6.0):
        msg = f"{sersic_index=}"
        mu = sersic_surfbright(radii, sersic_index=sersic_index)

        def func_integrate(r: float | np.ndarray[float]):
            return np.pi*r*sersic_surfbright(r=r, sersic_index=sersic_index)

        fluxes = sersic_integrated_flux(radii, sersic_index=sersic_index)
        for flux, r_min, r_max in zip(fluxes, radii_start, radii):
            flux_est, err, *_ = quad(func_integrate, r_min, r_max)
            np.testing.assert_almost_equal(
                flux, flux_est,
                decimal=max(err/flux, 1e-7),
                err_msg=f"integral for {sersic_index=} from {r_min=} to {r_max=}"
            )

        np.testing.assert_array_less(mu[1:], mu[:-1], err_msg=msg)
        np.testing.assert_array_less(mu, np.inf, err_msg=msg)

        mu_rescaled = sersic_surfbright(radii*1.05, sersic_index=sersic_index, r_eff=1.05)
        np.testing.assert_array_almost_equal(mu, mu_rescaled, decimal=8, err_msg=msg)
