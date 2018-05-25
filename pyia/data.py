# coding: utf-8
""" Data structures. """

# Standard library
from os import path

# Third-party
import astropy.coordinates as coord
from astropy.io import fits
from astropy.table import Table, Column
import astropy.units as u
import numpy as np
import pandas as pd

__all__ = ['GaiaData']


# This is from reading the data model
gaia_unit_map = {
    'ra': u.degree,
    'dec': u.degree,
    'parallax': u.milliarcsecond,
    'pmra': u.milliarcsecond / u.year,
    'pmdec': u.milliarcsecond / u.year,
    'radial_velocity': u.km / u.s,
    'ra_error': u.milliarcsecond,
    'dec_error': u.milliarcsecond,
    'parallax_error': u.milliarcsecond,
    'pmra_error': u.milliarcsecond / u.year,
    'pmdec_error': u.milliarcsecond / u.year,
    'radial_velocity_error': u.km / u.s,
    'astrometric_excess_noise': u.mas,
    'astrometric_weight_al': 1/u.mas**2,
    'astrometric_pseudo_colour': 1/u.micrometer,
    'astrometric_pseudo_colour_error': 1/u.micrometer,
    'astrometric_sigma5d_max': u.mas,
    'phot_g_mean_flux': u.photon/u.s,
    'phot_g_mean_flux_error': u.photon/u.s,
    'phot_g_mean_mag': u.mag,
    'phot_bp_mean_flux': u.photon/u.s,
    'phot_bp_mean_flux_error': u.photon/u.s,
    'phot_bp_mean_mag': u.mag,
    'phot_rp_mean_flux': u.photon/u.s,
    'phot_rp_mean_flux_error': u.photon/u.s,
    'phot_rp_mean_mag': u.mag,
    'bp_rp': u.mag,
    'bp_g': u.mag,
    'g_rp': u.mag,
    'rv_template_teff': u.K,
    'l': u.degree,
    'b': u.degree,
    'ecl_lon': u.degree,
    'ecl_lat': u.degree,
    'teff_val': u.K,
    'teff_percentile_lower': u.K,
    'teff_percentile_upper': u.K,
    'a_g_val': u.mag,
    'a_g_percentile_lower': u.mag,
    'a_g_percentile_upper': u.mag,
    'e_bp_min_rp_val': u.mag,
    'e_bp_min_rp_percentile_lower': u.mag,
    'e_bp_min_rp_percentile_upper': u.mag,
    'radius_val': u.Rsun,
    'radius_percentile_lower': u.Rsun,
    'radius_percentile_upper': u.Rsun,
    'lum_val': u.Lsun,
    'lum_percentile_lower': u.Lsun,
    'lum_percentile_upper': u.Lsun,
}


class GaiaData:
    """Class for loading and interacting with data from the Gaia mission. This
    should work with data from any data release, i.e., DR1 gaia_source or TGAS,
    or DR2 gaia_source.

    Parameters
    ----------
    data : `astropy.table.Table`, `pandas.DataFrame`, dict_like, str
        This must be pre-loaded data as any of the types listed above, or a
        string filename containing a table that is readable by
        `astropy.table.Table.read`.
    """

    def __init__(self, data, **kwargs):
        if not isinstance(data, Table):
            if isinstance(data, str):
                data = Table.read(data, **kwargs)

            else:
                # the dict-like object might have Quantity's, so we want to
                # preserve any units
                data = Table(data, **kwargs)

        # HACK and JFC: make sure table isn't masked
        if data.masked:
            cols = []
            for c in data.colnames:
                col = data[c]
                col.mask = None
                cols.append(Column(col))
            data = Table(cols)

        # Create a copy of the default unit map
        self.units = gaia_unit_map.copy()

        # Update the unit map with the table units
        for c in data.colnames:
            if data[c].unit is not None:
                self.units[c] = u.Unit(data[c].unit)

        # By this point, data should always be a DataFrame (for @smoh)
        self.data = data
        self._has_rv = 'radial_velocity' in self.data

        # For caching later
        self._cov = None
        self._cov_units = None
        self._coord = None

    @classmethod
    def from_query(cls, query_str, login_info=None):
        """
        Run the specified query and return a `GaiaData` instance with the
        returned data.

        This is meant only to be used for quick queries to the main Gaia science archive. For longer queries and more customized usage, use TAP access to any of the Gaia mirrors with, e.g., astroquery or pyvo.

        This requires ``astroquery`` to be installed.

        Parameters
        ----------
        query_str : str
            The string ADQL query to execute.
        login_info : dict, optional
            Username and password for the Gaia science archive as keys "user"
            and "password". If not specified, will use anonymous access, subject
            to the query limits.

        Returns
        -------
        gaiadata : `GaiaData`
            An instance of this object.

        """
        try:
            from astroquery.gaia import Gaia
        except ImportError:
            raise ImportError('Failed to import astroquery. To use the '
                              'from_query() classmethod, you must first'
                              ' install astroquery, e.g., with pip: '
                              '\n\tpip install astroquery')

        if login_info is not None:
            Gaia.login(**login_info)

        job = Gaia.launch_job_async(query_str)
        tbl = job.get_results()

        return cls(tbl)

    ##########################################################################
    # Python internal
    #
    def __getattr__(self, name):
        # to prevent recursion errors:
        # nedbatchelder.com/blog/201010/surprising_getattr_recursion.html
        if name in ['data', 'units']:
            raise AttributeError()

        if name in self.units:
            return np.asarray(self.data[name]) * self.units[name]

        else:
            return self.data[name]

    def __setattr__(self, name, val):

        if name in ['data', 'units']:
            # needs to be here to catch the first time we enter this func.
            super().__setattr__(name, val)

        elif name in self.units:
            if not hasattr(val, 'unit'):
                raise ValueError('To set data for column "{0}", you must '
                                 'provide a Quantity-like object (with units).'
                                 .format(name))
            self.data[name] = val
            self.units[name] = val.unit

        elif name in self.data.columns:
            self.data[name] = val

        else:
            super().__setattr__(name, val)

    def __dir__(self):
        return super().__dir__() + [str(k) for k in self.data.columns]

    def __getitem__(self, slc):
        if isinstance(slc, int):
            slc = slice(slc, slc+1)
        return self.__class__(self.data[slc])

    def __len__(self):
        return len(self.data)

    def __str__(self):
        names = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']
        if self._has_rv:
            names.append('radial_velocity')
        return str(self.data[names])

    def __repr__(self):
        return "<GaiaData: {0:d} rows>".format(len(self))

    ##########################################################################
    # Computed and convenience quantities
    #
    @property
    def pm(self):
        """2D proper motion. Has shape `(nrows, 2)`"""
        _u = self.pmra.unit
        return np.vstack((self.pmra.value, self.pmdec.to(_u).value)).T * _u

    @property
    def distance(self):
        """Assumes 1/parallax. Has shape `(nrows,)`"""
        return coord.Distance(parallax=self.parallax)

    @property
    def distmod(self):
        """Distance modulus, m-M = 5 * log10(dist / (10 pc))"""
        return self.distance.distmod

    @property
    def vtan(self):
        """
        Tangential velocity computed using the proper motion and inverse
        parallax as the distance. Has shape `(nrows, 2)`
        """
        d = self.distance
        vra = (self.pmra * d).to(u.km/u.s, u.dimensionless_angles()).value
        vdec = (self.pmdec * d).to(u.km/u.s, u.dimensionless_angles()).value
        return np.vstack((vra, vdec)).T * u.km/u.s

    def get_cov(self, RAM_threshold=1*u.gigabyte, units=None):
        """
        The Gaia data tables contain correlation coefficients and standard
        deviations for (ra, dec, parallax, pm_ra, pm_dec), but for most analyses
        we need covariance matrices. This converts the data provided by Gaia
        into covariance matrices.

        If a radial velocity exists, this also contains the radial velocity
        variance. If radial velocity doesn't exist, that diagonal element is set
        to inf.

        The default units of the covariance matrix are [degree, degree, mas,
        mas/yr, mas/yr, km/s], but this can be modified by passing in a
        dictionary with new units. For example, to change just the default ra,
        dec units for the covariance matrix, you can pass in::

            units=dict(ra=u.radian, dec=u.radian)

        Parameters
        ----------
        RAM_threshold : `astropy.units.Quantity`
            Raise an error if the expected covariance array is larger than the
            specified threshold. Set to ``None`` to disable this checking.
        """

        if RAM_threshold is not None:
            # Raise error if the user is going to blow up their RAM
            estimated_RAM = 6 * 6 * len(self) * 8*u.bit
            if estimated_RAM > RAM_threshold:
                raise RuntimeError('Estimated RAM usage for generating '
                                   'covariance matrices is larger than the '
                                   'specified threshold. Use the argument: '
                                   '`RAM_threshold=None` to disable this check')

        if units is None:
            units = dict()
        units.setdefault('ra', u.deg)
        units.setdefault('dec', u.deg)
        units.setdefault('parallax', u.mas)
        units.setdefault('pmra', u.mas/u.yr)
        units.setdefault('pmdec', u.mas/u.yr)
        units.setdefault('radial_velocity', u.km/u.s)

        # Use cached property if it exists
        if self._cov is not None and self._cov_units is not None:
            unit_check = [units[k] == self._cov_units[k] for k in units]
            if all(unit_check):
                return self._cov
        self._cov_units = units

        # The full returned matrix
        C = np.zeros((len(self), 6, 6))

        # We handle radial_velocity separately below - doesn't have correlation
        # coefficients with the astrometric parameters
        names = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']

        # pre-load the diagonal
        for i, name in enumerate(names):
            err = getattr(self, name + "_error")
            C[:, i, i] = err.to(units[name]).value ** 2

        if self._has_rv:
            name = 'radial_velocity'
            err = getattr(self, name + "_error")
            C[:, 5, 5] = err.to(units[name]).value ** 2
        else:
            C[:, 5, 5] = np.inf

        C[:, 5, 5][np.isnan(C[:, 5, 5])] = np.inf # missing values

        for i, name1 in enumerate(names):
            for j, name2 in enumerate(names):
                if j <= i:
                    continue

                corr = getattr(self, "{0}_{1}_corr".format(name1, name2))

                # We don't need to worry about units here because the diagonal
                # values have already been converted
                C[:, i, j] = corr * np.sqrt(C[:, i, i] * C[:, j, j])
                C[:, j, i] = C[:, i, j]

        self._cov = C
        return self._cov

    ##########################################################################
    # Astropy connections
    #
    @property
    def skycoord(self):
        """
        Return an `~astropy.coordinates.SkyCoord` object to represent
        all coordinates. Note: this requires Astropy v3.0 or higher!
        """
        return self.get_skycoord()

    def get_skycoord(self, distance=None, radial_velocity=None):
        """
        Return an `~astropy.coordinates.SkyCoord` object to represent
        all coordinates. Note: this requires Astropy v3.0 or higher!
        """
        # TODO: cache is disabled
        # if self._coord is None:

        kw = dict()
        if self._has_rv:
            kw['radial_velocity'] = self.radial_velocity

        if radial_velocity is not False and radial_velocity is not None:
            kw['radial_velocity'] = radial_velocity
        elif radial_velocity is False and 'radial_velocity' in kw:
            kw.pop('radial_velocity')

        if distance is None:
            kw['distance'] = self.distance
        elif distance is not False and distance is not None:
            kw['distance'] = distance

        _coord = coord.SkyCoord(ra=self.ra, dec=self.dec,
                                pm_ra_cosdec=self.pmra,
                                pm_dec=self.pmdec, **kw)

        return _coord
