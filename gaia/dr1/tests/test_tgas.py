# Standard library

# Third-party
import astropy.coordinates as coord
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
from astropy.tests.helper import quantity_allclose
import numpy as np

# Project
from ..tgas import TGASData, TGASStar

tgas_filename = get_pkg_data_filename('data/tiny-tgas.fits')

def test_init_dunder():
    """Test creating the object, and the dunder methods"""

    # Creation from a filename
    t1 = TGASData(tgas_filename)

    # Make sure creation from an existing table works:
    t2 = TGASData(fits.getdata(tgas_filename, 1))

    # __getattr__
    assert hasattr(t1, 'pm_ra_cosdec')
    t1.pm_ra_cosdec # make sure this works

    assert hasattr(t1, 'random_index')
    t1.random_index

    t1.data

    # __getitem__
    star = t1[0]
    assert isinstance(star, TGASStar)
    assert star.ra == t1.ra[0]

    t1_sub = t1[4:8]
    assert isinstance(t1_sub, TGASData)
    assert np.all(t1_sub.ra.value == t1.ra.value[4:8])

    idx = [1, 5, 99]
    t1_sub = t1[idx]
    assert isinstance(t1_sub, TGASData)
    assert np.all(t1_sub.ra.value == t1.ra.value[idx])

    # __dir__
    assert 'pm_ra_cosdec' in dir(t1)
    assert 'pmra' not in dir(t1)
    assert 'random_index' in dir(t1)

    # __len__
    assert len(t1) == 128

def test_convenience():
    tgas = TGASData(tgas_filename)

    # TODO: how to actually test below?

    # get_distance()
    d1 = tgas.get_distance()
    d2 = tgas.get_distance(lutz_kelker=True)
    assert not quantity_allclose(d1, d2)

    # get_vtan()
    v1 = tgas.get_vtan()
    v2 = tgas.get_vtan(lutz_kelker=True)
    assert not quantity_allclose(v1, v2)

    # get_coord()
    c1 = tgas.get_coord()
    c2 = tgas.get_coord(lutz_kelker=True)
    assert not quantity_allclose(c1.cartesian.xyz, c2.cartesian.xyz)

    # parallax_snr
    assert len(tgas.parallax_snr) == len(tgas)

    # distance
    assert isinstance(tgas.distance, coord.Distance)

def test_tgasstar():
    tgas = TGASData(tgas_filename)
    star = tgas[8]

    # make sure len() works:
    assert len(star) == 1

    # get_cov()
    cov = star.get_cov()
    assert cov.shape == (5,5)
    for i in range(5):
        for j in range(5):
            assert cov[i,j] == cov[j,i]

    # get_y()
    y = star.get_y()
    assert y.shape == (5,)

    # get_y_samples()
    y = star.get_y_samples()
    assert y.shape == (5,)

    y = star.get_y_samples(size=10)
    assert y.shape == (10,5)

    # get_coord_samples()
    c = star.get_coord_samples()
    assert c.shape == ()

    c = star.get_coord_samples(size=10)
    assert c.shape == (10,)
