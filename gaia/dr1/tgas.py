# Standard library
import os
from os import path
import sys

# Third-party
import astropy.coordinates as coord
from astropy.table import Table, Row
import astropy.units as u
from astropy.utils.compat.misc import override__dir__
import numpy as np

# Project
from .core import DataMeta, ComponentMapping

__all__ = ['TGASData', 'TGASStar']


class TGASData(metaclass=DataMeta):
    component_mappings = [
        ComponentMapping('ra', 'ra', u.degree),
        ComponentMapping('dec', 'dec', u.degree),
        ComponentMapping('parallax', 'parallax', u.mas),
        ComponentMapping('pm_ra_cosdec', 'pmra', u.mas/u.yr),
        ComponentMapping('pm_dec', 'pmdec', u.mas/u.yr),
        ComponentMapping('ra_error', 'ra_error', u.mas),
        ComponentMapping('dec_error', 'dec_error', u.mas),
        ComponentMapping('parallax_error', 'parallax_error', u.mas),
        ComponentMapping('pm_ra_cosdec_error', 'pmra_error', u.mas/u.yr),
        ComponentMapping('pm_dec_error', 'pmdec_error', u.mas/u.yr)
    ]

    def __init__(self, filename_or_data):

        if isinstance(filename_or_data, str):
            self.data = Table.read(filename_or_data, hdu=1)

        else:
            self.data = Table(filename_or_data)

    # --------------------------------------------------------------------------
    # Python special methods:
    #
    def __getattr__(self, name):
        # to prevent recursion errors:
        # http://nedbatchelder.com/blog/201010/surprising_getattr_recursion.html
        if name == 'data':
            raise AttributeError()

        if name in self._obj_to_tbl:
            tbl_name = self._obj_to_tbl[name]
            unit = self._obj_to_unit[name]
            return self.data[tbl_name] * unit

        elif name in self.data.colnames:
            return self.data[name]

        else:
            raise AttributeError("Object {} has no attribute '{}' and source "
                                 "data table has no column with that name."
                                 .format(self, name))

    def __getitem__(self, slc):
        sliced = self.data[slc]

        if isinstance(sliced, Row): # this is only one row
            return TGASStar(sliced)

        else: # multiple rows
            return TGASData(sliced)

    @override__dir__
    def __dir__(self):
        """Override the builtin `dir` behavior to include representation names.
        """
        colnames = []
        for name in self.data.colnames:
            if (name in self._tbl_to_obj and
                    self._tbl_to_obj[name] not in self._extra_dir_names):
                colnames.append(name)

            elif name not in self._tbl_to_obj:
                colnames.append(name)

        return self._extra_dir_names + colnames

    def __len__(self):
        return len(self.data)

    # --------------------------------------------------------------------------
    # Other convenience methods:
    #
    def get_distance(self, lutz_kelker=False):
        """Return the distance with or without the Lutz-Kelker correction.
        """

        if lutz_kelker:
            snr = self.parallax / self.parallax_error
            tmp = self.parallax * (0.5 + 0.5*np.sqrt(1 - 16/snr**2))

        else:
            tmp = self.parallax

        return coord.Distance(tmp.to(u.pc, u.parallax()))

    def get_vtan(self, lutz_kelker=False):
        """Return the tangential velocity computed using the proper motion and
        distance.
        """
        d = self.get_distance(lutz_kelker=lutz_kelker)
        vra = (self.pm_ra_cosdec * d).to(u.km/u.s,
                                         u.dimensionless_angles()).value
        vdec = (self.pm_dec * d).to(u.km/u.s, u.dimensionless_angles()).value
        return np.vstack((vra, vdec)) * u.km/u.s

    def get_coord(self, lutz_kelker=False):
        """
        Return an `astropy.coordinates` object to represent all coordinates.
        """
        return coord.ICRS(ra=self.ra, dec=self.dec,
                          distance=self.get_distance(lutz_kelker=lutz_kelker),
                          pm_ra_cosdec=self.pm_ra_cosdec,
                          pm_dec=self.pm_dec)

    # --------------------------------------------------------------------------
    # Custom attributes:
    #
    @property
    def parallax_snr(self):
        """Parallax signal-to-noise ratio."""
        return self.parallax / self.parallax_error

    @property
    def distance(self):
        """Distance, computed as 1/parallax."""
        return self.get_distance(lutz_kelker=False)

    # --------------------------------------------------------------------------
    # Methods that call teh internets:
    #
    @classmethod
    def download(cls):
        """Download and cache the TGAS catalog.

        Parameters
        ----------
        """
        pass

    # TODO: some wrappers around astroquery calls to get other catalog
    # information and join it in to the data table


class TGASStar(TGASData):

    def __init__(self, row):
        self.data = row
        self._cov = None # for caching
        self._Cinv = None # for caching

    def __len__(self):
        return 1

    def get_cov(self):
        """Retrieve the 5x5 covariance matrix for all astrometric data
        components in the units of the components.
        """

        tbl_names = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']
        if self._cov is None:
            C = np.zeros((5, 5))

            # pre-load the diagonal
            for i,name in enumerate(tbl_names):
                err_name = "{0}_error".format(name)
                fac = self._err_unit_scale_factor[err_name]
                C[i,i] = (self.data[err_name] * fac) ** 2

            for i,name1 in enumerate(tbl_names):
                for j,name2 in enumerate(tbl_names):
                    if j <= i: # we make it symmetric below
                        continue

                    full_name = "{0}_{1}_corr".format(name1, name2)
                    if full_name not in self.data.colnames:
                        continue # skip if no correlations exist

                    # Don't need to deal with units because the diagonals have
                    # already been converted to the desired units
                    C[i,j] = self.data[full_name] * np.sqrt(C[i,i]*C[j,j])
                    C[j,i] = self.data[full_name] * np.sqrt(C[i,i]*C[j,j])

            self._cov = C

        return self._cov

    def get_y(self):
        y = np.array([self.ra.value, self.dec.value, self.parallax.value,
                      self.pm_ra_cosdec.value, self.pm_dec.value])
        return y

    def get_y_samples(self, size=1):
        y = self.get_y()
        Cov = self.get_cov()
        samples = np.random.multivariate_normal(y, Cov, size=size)
        return samples

    def get_coord_samples(self, size=1):
        """Create an `astropy.coordinates` object with samples from the error
        distribution over astrometric data.
        """
        samples = self.get_y_samples(self, size=size)
        return coord.ICRS(ra=samples[:,0] * u.deg,
                          dec=samples[:,1] * u.deg,
                          distance=1000./samples[:,2] * u.pc,
                          pm_ra_cosdec=samples[:,3] * u.mas/u.yr,
                          pm_dec=samples[:,4] * u.mas/u.yr)
