import lsst.gauss2d.fit as g2f

import pytest


def test_sersicmix():
    available = g2f.sersic_mix_descriptions()

    assert len(available) > 0

    for order, versions in available.items():
        assert len(versions) > 0
        for version, description in versions.items():
            assert len(version) > 0
            assert len(description) > 0

            knots = g2f.sersic_mix_knots(order, version)
            assert len(knots) > 0