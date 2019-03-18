# Standard library
from os import path
import numpy as np

# Third-party
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename
import pytest

# Project
from ..data import GaiaData, gaia_unit_map

@pytest.mark.parametrize('filename,fmt', [
    (get_pkg_data_filename('data/gdr2_sm.fits'), 'fits'),
    (get_pkg_data_filename('data/gdr2_sm.vot'), 'votable')
])
def test_init(filename, fmt):
    """Test initializing the GaiaData class with different options"""

    if 'fits' in path.splitext(filename)[1]:
        gd = GaiaData(fits.getdata(filename, 1))
        assert isinstance(gd.data, Table)

    print(filename)
    gd = GaiaData(Table.read(filename, format=fmt))
    gd = GaiaData(filename, format=fmt)


def test_slicing_getattr():
    filename = get_pkg_data_filename('data/gdr2_sm.fits')
    gd = GaiaData(filename)

    for k, unit in gaia_unit_map.items():
        if k in gd.data.columns:
            # have to special case the flux columns:
            if unit == u.ph/u.s:
                continue
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
    gd.a_g_val


def test_str_repr():
    filename = get_pkg_data_filename('data/gdr2_sm.fits')
    gd = GaiaData(Table.read(filename))
    assert 'GaiaData' in repr(gd)
    assert '100 rows' in repr(gd)


def test_computed_quantities():
    filename = get_pkg_data_filename('data/gdr2_sm.fits')
    gd = GaiaData(Table.read(filename))

    assert gd.pm.unit == u.mas/u.yr
    assert gd.pm.shape == (len(gd), 2)

    assert isinstance(gd.distance, coord.Distance)
    assert gd.distance.shape == (len(gd), )

    assert gd.distmod.unit == u.mag

    assert gd.vtan.unit == u.km/u.s
    assert gd.vtan.shape == (len(gd), 2)


def test_cov():
    filename = get_pkg_data_filename('data/gdr2_sm.fits')
    gd = GaiaData(Table.read(filename))

    C = gd.get_cov()
    assert C.shape == (len(gd), 6, 6)

    with pytest.raises(RuntimeError):
        gd.get_cov(RAM_threshold=1*u.kilobyte) # force failure


def test_skycoord():
    filename = get_pkg_data_filename('data/gdr2_sm.fits')
    gd = GaiaData(Table.read(filename))

    c = gd.skycoord
    assert len(c) == len(gd)


def test_setattr():
    filename = get_pkg_data_filename('data/gdr2_sm.fits')
    tbl = Table.read(filename)
    tbl['arr_column'] = np.arange(len(tbl)) * 10.
    gd = GaiaData(tbl)

    # Setting a quantity column
    with pytest.raises(ValueError):
        gd.parallax = np.arange(len(gd))

    new_vals = np.arange(len(gd)) * u.microarcsecond
    gd.parallax = new_vals
    assert np.all(np.asarray(gd.data['parallax']) == new_vals.value)
    assert gd.parallax.unit == u.microarcsecond

    # Setting an array column
    new_vals = np.arange(len(gd))
    gd.arr_column = new_vals
    assert np.all(gd.data['arr_column'] == new_vals)


@pytest.mark.remote_data
@pytest.mark.xfail # Gaia archive is down!
def test_from_query():
    q = '''SELECT TOP 10 * FROM gaiadr1.tgas_source'''
    gd = GaiaData.from_query(q)

    assert len(gd) == 10


def test_get_samples():
    filename = get_pkg_data_filename('data/gdr2_sm.fits')
    gd = GaiaData(Table.read(filename))

    g_samples = gd.get_error_samples(size=16)
    assert g_samples.ra.shape == (len(gd), 16)
    assert g_samples.dec.shape == (len(gd), 16)
    assert g_samples.parallax.shape == (len(gd), 16)
    assert g_samples.pmra.shape == (len(gd), 16)
    assert g_samples.pmdec.shape == (len(gd), 16)
    assert g_samples.radial_velocity.shape == (len(gd), 16)

    c = g_samples.get_skycoord(distance=False)
    assert c.shape == (100, 16)
