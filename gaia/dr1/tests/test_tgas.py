# Standard library

# Third-party
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import numpy as np

# Project
from ..tgas import TGASData, TGASStar

def test_init_dunder():
    """Test creating the object, and the dunder methods"""
    tgas_filename = get_pkg_data_filename('data/tiny-tgas.fits')

    t1 = TGASData(tgas_filename)
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
