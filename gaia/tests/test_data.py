# Standard library
from os import path
import numpy as np

# Third-party
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename
import pandas as pd
import pytest

# Project
from ..data import GaiaData, gaia_unit_map

@pytest.mark.parametrize('filename', [
    get_pkg_data_filename('data/gdr2mock.fits'),
    get_pkg_data_filename('data/gdr2mock_unitless.fits'),
    get_pkg_data_filename('data/gdr2mock.votable')
])
def test_init(filename):
    """Test initializing the GaiaData class with different options"""

    if 'fits' in path.splitext(filename)[1]:
        gd = GaiaData(fits.getdata(filename, 1))
        assert isinstance(gd.data, pd.DataFrame)

    gd = GaiaData(Table.read(filename))
    gd = GaiaData(Table.read(filename).to_pandas())


def test_slicing_getattr():
    filename = get_pkg_data_filename('data/gdr2mock.fits')
    gd = GaiaData(Table.read(filename))

    for k, unit in gaia_unit_map.items():
        assert getattr(gd, k).unit == unit

    size = len(gd)
    gd_slc = gd[:4]
    assert len(gd_slc) == 4

    gd_one = gd[4]
    assert len(gd_one) == 1

    assert 'parallax' in dir(gd)
    assert 'source_id' in dir(gd)
    gd.source_id
    gd.parallax
    gd.ra_dec_corr


def test_str_repr():
    filename = get_pkg_data_filename('data/gdr2mock.fits')
    gd = GaiaData(Table.read(filename))
    assert 'GaiaData' in repr(gd)
    assert '100 rows' in repr(gd)


def test_computed_quantities():
    filename = get_pkg_data_filename('data/gdr2mock.fits')
    gd = GaiaData(Table.read(filename))

    assert gd.pm.unit == u.mas/u.yr
    assert gd.pm.shape == (len(gd), 2)

    assert isinstance(gd.distance, coord.Distance)
    assert gd.distance.shape == (len(gd), )

    assert gd.distmod.unit == u.mag

    assert gd.vtan.unit == u.km/u.s
    assert gd.vtan.shape == (len(gd), 2)


def test_cov():
    filename = get_pkg_data_filename('data/gdr2mock.fits')
    gd = GaiaData(Table.read(filename))

    C = gd.get_cov()
    assert C.shape == (len(gd), 6, 6)

    with pytest.raises(RuntimeError):
        gd.get_cov(RAM_threshold=1*u.kilobyte) # force failure


def test_skycoord():
    filename = get_pkg_data_filename('data/gdr2mock.fits')
    gd = GaiaData(Table.read(filename))

    c = gd.skycoord
    assert len(c) == len(gd)
