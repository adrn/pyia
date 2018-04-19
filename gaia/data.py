# coding: utf-8
""" Data structures. """

# Standard library
import numpy as np

# Third-party
import astropy.coordinates as coord
from astropy.table import Table
import astropy.units as u
import pandas as pd

# Project
from ..utils import cached_property

__all__ = ['GaiaData']


# This is from reading the data model
gaia_unit_map = {
    'ra': u.degree,
    'dec': u.degree,
    'parallax': u.milliarcsecond,
    'pmra': u.milliarcsecond / u.year,
    'pmdec': u.milliarcsecond / u.year,
    'radial_velocity': u.km / u.s,
    'ra_error': u.degree,
    'dec_error': u.degree,
    'parallax_error': u.milliarcsecond,
    'pmra_error': u.milliarcsecond / u.year,
    'pmdec_error': u.milliarcsecond / u.year,
    'radial_velocity_error': u.km / u.s
}


class GaiaData:
    """Class for loading and interacting with data from the Gaia mission. This
    should work with data from any data release, i.e., DR1 gaia_source or TGAS,
    or DR2 gaia_source.

    Parameters
    ----------
    data : `astropy.table.Table`, `pandas.DataFrame`, dict_like
        This must be pre-loaded data as any of the types listed above.
    """

    def __init__(self, data):
        if not isinstance(data, Table) or not isinstance(data, pd.DataFrame):
            # the dict-like object might have Quantity's, so we want to
            # preserve any units
            data = Table(data)

        # Create a copy of the default unit map
        self.units = gaia_unit_map.copy()

        if isinstance(data, Table):
            # Modify unit dict if the input object has custom units:
            for name in data.columns:
                if hasattr(data[name], 'unit'):
                    self.units[name] = data[name].unit
            data = data.to_pandas()

        # By this point, data should always be a DataFrame (for @smoh)
        self.data = data

    @classmethod
    def from_query(cls, query_str, ):
        pass


    ##########################################################################
    # Python internal
    #
    def __getattr__(self, name):
        # to prevent recursion errors:
        # nedbatchelder.com/blog/201010/surprising_getattr_recursion.html
        if name == 'data':
            raise AttributeError()

        if name in GaiaData.units:
            return self.data[name].values * self.units[name]

        else:
            return self.data[name]

    def __getitem__(self, slc):
        sliced = self.data[slc]

        if sliced.ndim == 0: # only one row after slice
            return GaiaSource(sliced)

        else: # multiple rows
            return GaiaData(sliced)

    def __len__(self):
        return len(self.data)

    ##########################################################################
    # Computed quantities
    #
    def get_vtan(self, lutz_kelker=True):
        """
        Return the tangential velocity computed using the proper motion
        and distance.
        """
        d = self.get_distance(lutz_kelker=lutz_kelker)
        vra = (self.pmra * d).to(u.km/u.s, u.dimensionless_angles()).value
        vdec = (self.pmdec * d).to(u.km/u.s, u.dimensionless_angles()).value
        return np.vstack((vra, vdec)).T * u.km/u.s

    ##########################################################################
    # Astropy connections
    #
    @cached_property
    def coords(self):
        """
        Return an `~astropy.coordinates.SkyCoord` object to represent
        all coordinates.
        """
        return coord.SkyCoord(ra=self.ra, dec=self.dec,
                              distance=self.get_distance(lutz_kelker=lutz_kelker))

    @property
    def distmod(self):
        """Distance modulus, m-M = 5 * log10(dist / (10 pc))"""
        return coord.Distance(parallax=self.parallax).distmod

    @property
    def parallax_snr(self):
        return self.parallax / self.parallax_error



class TGASStar(TGASData):
    def __init__(self, row, rv=None, rv_err=None):
        self._data = row
        self._cov = None # for caching
        self._Cinv = None # for caching

        # radial velocities
        if rv is not None:
            if rv_err is None:
                raise ValueError("If radial velocity is provided, you must also "
                                 "provide an error.")
            self._rv = rv.to(u.km/u.s).value
            self._rv_err = rv_err.to(u.km/u.s).value
        else:
            self._rv = 0. # need to do this so y can have float type
            self._rv_err = None

    def __len__(self):
        return 1

    def __getitem__(self, slc):
        object.__getitem__(self, slc)

    def __str__(self):
        infostr = '\n'.join([
            # 'index    = %i' %(i),
            'ra       = %s' % (self.ra),
            'dec      = %s' % (self.dec),
            'parallax = %s (snr = %.1f)' % (self.parallax, self.parallax_snr),
            'pmra     = %s (snr = %.1f)' % (self.pmra, self.pmra/self.pmra_error),
            'pmdec    = %s (snr = %.1f)' % (self.pmdec, self.pmdec/self.pmdec_error),
            'dist vra vdec = %s %s' % (self.get_distance(), self.get_vtan()),
        ])
        return infostr

    def get_cov(self):
        """
        The Gaia TGAS data table contains correlation coefficients and standard
        deviations for (ra, dec, parallax, pm_ra, pm_dec), but for most analysis
        we need covariance matrices. This converts the Gaia table into covariance
        matrix. If a radial velocity was specified on creation, this also contains
        the radial velocity variance. The base units are:
        [deg, deg, mas, mas/yr, mas/yr, km/s]
        """

        if self._cov is not None:
            return self._cov

        names = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']

        C = np.zeros((6,6))

        # pre-load the diagonal
        for i,name in enumerate(names):
            full_name = "{}_error".format(name)
            C[i,i] = self._data[full_name]**2

        for i,name1 in enumerate(names):
            for j,name2 in enumerate(names):
                if j <= i:
                    continue
                full_name = "{}_{}_corr".format(name1, name2)
                C[i,j] = self._data[full_name] * np.sqrt(C[i,i]*C[j,j])
                C[j,i] = self._data[full_name] * np.sqrt(C[i,i]*C[j,j])

        if self._rv_err is not None:
            C[5,5] = self._rv_err**2

        self._cov = C
        return self._cov



def angdist_to(ra1, dec1, ra2, dec2):
    """Return angular distance in degrees"""
    ra1, dec1 = np.deg2rad(ra1), np.deg2rad(dec1)
    ra2, dec2 = np.deg2rad(ra2), np.deg2rad(dec2)
    return np.rad2deg(np.arccos(np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)))
