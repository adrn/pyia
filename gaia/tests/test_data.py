# Standard library
import numpy as np

# Third-party
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u
from astropy.utils.data import get_pkg_data_fileobj

# Project
from ..core import GaiaData


def test_init():
    """Test initializing the GaiaData class with different options"""
    fn = get_pkg_data_fileobj('data/gdr2mock.fits')
