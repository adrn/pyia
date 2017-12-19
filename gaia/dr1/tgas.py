# Standard library
import os
from os import path
import sys

# Third-party
import astropy.coordinates as coord
from astropy.table import Table, Row
import astropy.units as u
import numpy as np

# Project
from .core import ComponentMapping, BaseGaiaData

__all__ = ['TGASData', 'TGASStar']


class TGASData(BaseGaiaData):
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
    def __getitem__(self, slc):
        sliced = self.data[slc]

        if isinstance(sliced, Row): # this is only one row
            return TGASStar(sliced)

        else: # multiple rows
            return TGASData(sliced)

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

        For example, in the source TGAS file, the error in RA is given in mas,
        but the value of RA is given in degrees. Here, the covariance matrix
        elements that include the RA error will have "degree" units, not mas.

        Returns
        -------
        cov : `numpy.ndarray`
            A 5 by 5 array containing the covariance matrix for all astrometric
            components for this star.
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
        """Get a single data vector containing RA, Dec, parallax, RA proper
        motion, Dec proper motion. The RA proper motion already includes the
        ``cos(dec)`` term.

        Returns
        -------
        y : `numpy.ndarray`
            A length-5 array of astrometric data.
        """
        y = np.array([self.ra.value, self.dec.value, self.parallax.value,
                      self.pm_ra_cosdec.value, self.pm_dec.value])
        return y

    def get_y_samples(self, size=None):
        """Get samples from the error distribution over the data vector (RA,
        Dec, parallax, RA proper motion, Dec proper motion). The RA proper
        motion already includes the ``cos(dec)`` term.

        Parameters
        ----------
        size : int, optional
            The number of samples to generate.

        Returns
        -------
        y_samples : `numpy.ndarray`
            A shape ``(size, 5)`` array of samples from the astrometric data
            error distribution.
        """
        y = self.get_y()
        Cov = self.get_cov()
        samples = np.random.multivariate_normal(y, Cov, size=size)
        return samples

    def get_coord_samples(self, size=None):
        """Create an `astropy.coordinates` object with samples from the error
        distribution over astrometric data.

        Parameters
        ----------
        size : int, optional
            The number of samples to generate.

        Returns
        -------
        coord_samples : `~astropy.coordinates.ICRS`
            A single array-valued `astropy.coordinates` object containing the
            samples from the error distribution.

        Examples
        --------
        TODO: combining with radial velocity data

        """
        samples = self.get_y_samples(size=size).T
        return coord.ICRS(ra=samples[0] * u.deg,
                          dec=samples[1] * u.deg,
                          distance=1000./samples[2] * u.pc,
                          pm_ra_cosdec=samples[3] * u.mas/u.yr,
                          pm_dec=samples[4] * u.mas/u.yr)
