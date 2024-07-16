import lsst.gauss2d.fit as g2f

import pytest


def test_channel():
    assert (len(g2f.Channel.all) == 1)

    with pytest.raises(ValueError):
        g2f.Channel("None")
    x = g2f.Channel("x")
    assert (len(g2f.Channel.all) == 2)
    with pytest.raises(ValueError):
        g2f.Channel.erase("None")
    with pytest.raises(ValueError):
        g2f.Channel("x")
    del x
    g2f.Channel.erase("x")
    assert (len(g2f.Channel.all) == 1)
    assert (len(g2f.Channel("x").all) == 2)
    assert (g2f.Channel.find("y") is None)
    assert (len(g2f.Channel("y").all) == 3)
