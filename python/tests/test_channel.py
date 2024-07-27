import lsst.gauss2d.fit as g2f

import pytest


def test_channel():
    # Previous tests may have added channels
    # Ideally we'd clean up but it is more difficult to ensure fixtures remove
    # all channels and uses thereof than to just accept that some may exist
    channels_all = g2f.Channel.all
    n_channels = len(channels_all)
    assert n_channels >= 1

    with pytest.raises(ValueError):
        # The None channel is special and accessed as g2f.Channel.NONE
        g2f.Channel("None")

    name_new = "A long name that is very unlikely to have already been taken"
    channel_new = g2f.Channel(name_new)
    assert len(g2f.Channel.all) == (n_channels + 1)
    with pytest.raises(ValueError):
        g2f.Channel.erase("None")
    with pytest.raises(ValueError):
        g2f.Channel(name_new)

    # deleting a Channel object should be fine
    del channel_new
    # ...only calling erase removes it from the registry, though
    g2f.Channel.erase(name_new)
    assert g2f.Channel.find(name_new) is None
    assert len(g2f.Channel.all) == n_channels
    assert len(g2f.Channel(name_new).all) == (n_channels + 1)
    assert len(g2f.Channel("Another unusually long name").all) == (n_channels + 2)
