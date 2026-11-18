# This file is part of gauss2d.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from typing import Iterable

import lsst.gauss2d.fit
import numpy as np
from scipy.special import gammainc
from scipy.stats import gamma as gamma_func

__all__ = [
    "gaussian_integrated_flux", "gaussian_mixture_integrated_flux", "sersic_integrated_flux",
    "sersic_surfbright",
]


def gaussian_integrated_flux(r: float | np.ndarray[float], sigma: float = 1.0):
    """Return the fraction of the integrated flux within given radii for a
    Gaussian profile.

    Parameters
    ----------
    r
        Radius or radii to compute the integrated flux within.
    sersic_index
        Sersic index of the profile.
    r_eff
        Effective radius of the profile.

    Returns
    -------
    flux_frac
        The fraction of integrated flux within the values of r.
    """
    if not (sigma > 0):
        raise ValueError(f"{sigma=} must be >0")
    return 1 - np.exp(-0.5*(r/sigma)**2)


def gaussian_mixture_integrated_flux(
    r: float | np.ndarray[float],
    integralsizes: Iterable[lsst.gauss2d.fit.IntegralSize],
):
    """Return the integrated flux within given radii for a Gaussian mixture.

    Parameters
    ----------
    r
        Radius or radii to compute the integrated flux within.
    sersic_index
        Sersic index of the profile.
    r_eff
        Effective radius of the profile.

    Returns
    -------
    fluxes
        The integrated flux within the values of r.
    """
    result = None
    for integralsize in integralsizes:
        flux_fracs = integralsize.integral*gaussian_integrated_flux(r=r, sigma=integralsize.sigma)
        if result is None:
            result = flux_fracs
        else:
            result += flux_fracs
    return result


def sersic_integrated_flux(r: float | np.ndarray[float], sersic_index: float, r_eff: float = 1.0):
    """Return the fraction of the integrated flux within given radii for a
    Sersic profile.

    Parameters
    ----------
    r
        Radius or radii to compute the integrated flux within.
    sersic_index
        Sersic index of the profile.
    r_eff
        Effective radius of the profile.

    Returns
    -------
    flux_frac
        The fraction of integrated flux within the values of r.
    """
    if not (sersic_index > 0) or not (r_eff > 0):
        raise ValueError(f"{sersic_index=} and {r_eff=} must be >0")
    two_ns = 2*sersic_index
    bn = gamma_func.ppf(0.5, two_ns)
    x = bn*(r/r_eff)**(1/sersic_index)
    flux_frac = gammainc(two_ns, x)
    return flux_frac


def sersic_surfbright(r: float | np.ndarray[float], sersic_index: float, r_eff: float = 1.0):
    """Return the surface brightness of a unit Sersic profile.

    Parameters
    ----------
    r
        Radius or radii to compute the surface brightness for.
    sersic_index
        Sersic index of the profile.
    r_eff
        Effective radius of the profile.

    Returns
    -------
    mu
        The surface brightness at the values of r.
    """
    if not (sersic_index > 0) or not (r_eff > 0):
        raise ValueError(f"{sersic_index=} and {r_eff=} must be >0")
    bn = gamma_func.ppf(0.5, 2*sersic_index)
    ir = np.exp(-bn*((r/r_eff)**(1./sersic_index) - 1.))
    return ir
