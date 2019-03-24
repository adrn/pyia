# coding: utf-8
""" Data structures. """

# Third-party
import astropy.coordinates as coord
from astropy.table import Table, Column
from astropy.time import Time
import astropy.units as u
import numpy as np

from .extinction import get_ext

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
    'ref_epoch': u.year
}

DEFAULT_REF_EPOCH = Time(2015.5, format='decimalyear')


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
            data = Table(cols, copy=False)

        # Create a copy of the default unit map
        self.units = gaia_unit_map.copy()

        # Store the source table
        self.data = data

        # Update the unit map with the table units
        self._invalid_units = dict()
        for c in data.colnames:
            if data[c].unit is not None:
                try:
                    self.units[c] = u.Unit(str(data[c].unit))
                except ValueError:
                    self._invalid_units[c] = data[c].unit

        self._has_rv = 'radial_velocity' in self.data.colnames

        # For caching later
        self._cache = dict()

    @classmethod
    def from_query(cls, query_str, login_info=None):
        """
        Run the specified query and return a `GaiaData` instance with the
        returned data.

        This is meant only to be used for quick queries to the main Gaia science
        archive. For longer queries and more customized usage, use TAP access to
        any of the Gaia mirrors with, e.g., astroquery or pyvo.

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

    @u.quantity_input(min_parallax=u.mas, equivalencies=u.parallax())
    def get_distance(self, min_parallax=None, parallax_fill_value=np.nan,
                     allow_negative=False):
        """Compute distance from parallax (by inverting the parallax) using
        `~astropy.coordinates.Distance`.

        Parameters
        ----------
        min_parallax : `~astropy.units.Quantity` (optional)
            If `min_parallax` specified, the parallaxes are clipped to this
            values (and it is also used to replace NaNs).
        allow_negative : bool (optional)
            This is passed through to `~astropy.coordinates.Distance`.

        Returns
        -------
        dist : `~astropy.coordinates.Distance`
            A ``Distance`` object with the data.
        """

        plx = self.parallax.copy()

        if np.isnan(parallax_fill_value):
            parallax_fill_value = parallax_fill_value * u.mas

        if min_parallax is not None:
            clipped = plx < min_parallax
            clipped |= ~np.isfinite(plx)
            plx[clipped] = parallax_fill_value

        return coord.Distance(parallax=plx, allow_negative=allow_negative)

    @property
    def distance(self):
        """Assumes 1/parallax. Has shape `(nrows,)`.

        This attribute will raise an error when there are negative or zero
        parallax values. For more flexible retrieval of distance values and
        auto-filling bad values, use the .get_distance() method."""
        return self.get_distance()

    def get_radial_velocity(self, fill_value=None):
        """Return radial velocity but with invalid values filled with the
        specified fill value.

        Parameters
        ----------
        fill_value : `~astropy.units.Quantity` (optional)
            If not ``None``, fill any invalid values with the specified value.
        """
        rv = self.radial_velocity.copy()
        rv[~np.isfinite(rv)] = fill_value
        return rv

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

        if 'cov' in self._cache:
            if units == self._cache['cov_units']:
                return self._cache['cov']

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

        # The full returned matrix
        C = np.zeros((len(self), 6, 6))

        # We handle radial_velocity separately below - doesn't have correlation
        # coefficients with the astrometric parameters
        names = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']

        # pre-load the diagonal
        for i, name in enumerate(names):
            if name + "_error" in self.data.colnames:
                err = getattr(self, name + "_error")
                C[:, i, i] = err.to(units[name]).value ** 2
            else:
                C[:, i, i] = np.nan

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

                if "{0}_{1}_corr".format(name1, name2) in self.data.colnames:
                    corr = getattr(self, "{0}_{1}_corr".format(name1, name2))
                else:
                    corr = np.nan

                # We don't need to worry about units here because the diagonal
                # values have already been converted
                C[:, i, j] = corr * np.sqrt(C[:, i, i] * C[:, j, j])
                C[:, j, i] = C[:, i, j]

        self._cache['cov'] = C
        self._cache['cov_units'] = units

        return self._cache['cov']

    def get_ebv(self, dustmaps_cls=None):
        """Compute the E(B-V) reddening at this location

        This requires the `dustmaps <http://dustmaps.readthedocs.io>`_ package
        to run!

        Parameters
        ----------
        dustmaps_cls : ``dustmaps`` query class
            By default, ``SFDQuery``.
        """
        if dustmaps_cls is None:
            from dustmaps.sfd import SFDQuery
            dustmaps_cls = SFDQuery

        c = self.get_skycoord(distance=False)
        return dustmaps_cls().query(c)

    def get_ext(self, dustmaps_cls=None):
        """Compute the E(B-V) reddening at this location

        This requires the `dustmaps <http://dustmaps.readthedocs.io>`_ package
        to run!

        Parameters
        ----------
        dustmaps_cls : ``dustmaps`` query class
            By default, ``SFDQuery``.

        Returns
        -------
        A_G
        A_BP
        A_RP
        """
        if 'A_G' not in self._cache:
            EBV = self.get_ebv(dustmaps_cls=dustmaps_cls)
            A_G, A_B, A_R = get_ext(self.phot_g_mean_mag.value,
                                    self.phot_bp_mean_mag.value,
                                    self.phot_rp_mean_mag.value,
                                    EBV)

            self._cache['A_G'] = A_G * u.mag
            self._cache['A_B'] = A_B * u.mag
            self._cache['A_R'] = A_R * u.mag

        return (self._cache['A_G'],
                self._cache['A_B'],
                self._cache['A_R'])

    def get_G0(self):
        """Return the extinction-corrected G-band magnitude."""
        A, _, _ = self.get_ext()
        return self.phot_g_mean_mag - A

    def get_BP0(self):
        """Return the extinction-corrected G_BP magnitude."""
        _, A, _ = self.get_ext()
        return self.phot_bp_mean_mag - A

    def get_RP0(self):
        """Return the extinction-corrected G_RP magnitude."""
        _, _, A = self.get_ext()
        return self.phot_rp_mean_mag - A

    ##########################################################################
    # Astropy connections
    #
    @property
    def skycoord(self):
        """
        Return an `~astropy.coordinates.SkyCoord` object to represent
        all coordinates. Note: this requires Astropy v3.0 or higher!

        Use the ``get_skycoord()`` method for more flexible access.
        """
        return self.get_skycoord()

    def get_skycoord(self, distance=None, radial_velocity=None,
                     ref_epoch=DEFAULT_REF_EPOCH):
        """
        Return an `~astropy.coordinates.SkyCoord` object to represent
        all coordinates. Note: this requires Astropy v3.0 or higher!

        `ref_epoch` is used to set the `obstime` attribute on the coordinate
        objects.  This is often included in the data release tables, but
        `ref_epoch` here is used if it's not.

        Parameters
        ----------
        distance : `~astropy.coordinate.Distance`, `~astropy.units.Quantity`, ``False`` (optional)
            If ``None``, this inverts the parallax to get the distance from the
            Gaia data. If ``False``, distance information is ignored. If an
            astropy ``Quantity`` or ``Distance`` object, it sets the distance
            values of the output ``SkyCoord`` to whatever is passed in.
        radial_velocity : `~astropy.units.Quantity` (optional)
            If ``None``, this uses radial velocity data from the input Gaia
            table. If an astropy ``Quantity`` object, it sets the radial
            velocity values of the output ``SkyCoord`` to whatever is passed in.
        ref_epoch : `~astropy.time.Time`, float (optional)
            The reference epoch of the data. If not specified, this will try to
            read it from the input Gaia data table. If not provided, this will
            be set to whatever the most recent data release is, so, **beware**!

        Returns
        -------
        c : `~astropy.coordinates.SkyCoord`
            The coordinate object constructed from the input Gaia data.
        """
        _coord_opts = (distance, radial_velocity)
        if 'coord' in self._cache:
            try:
                _check = self._cache['coord_opts'] == _coord_opts
            except ValueError: # array passed in for distance or radial_velocity
                _check = False

            if _check:
                return self._cache['coord']

        kw = dict()
        if self._has_rv:
            kw['radial_velocity'] = self.radial_velocity

        # Reference epoch
        if 'ref_epoch' in self.data.colnames:
            obstime = Time(self.ref_epoch.value, format='decimalyear')
        else:
            obstime = Time(ref_epoch, format='decimalyear')

        kw['obstime'] = obstime

        if radial_velocity is not False and radial_velocity is not None:
            kw['radial_velocity'] = radial_velocity
        elif radial_velocity is False and 'radial_velocity' in kw:
            kw.pop('radial_velocity')

        if distance is None:
            kw['distance'] = self.distance
        elif distance is not False and distance is not None:
            kw['distance'] = distance

        self._cache['coord'] = coord.SkyCoord(ra=self.ra, dec=self.dec,
                                              pm_ra_cosdec=self.pmra,
                                              pm_dec=self.pmdec, **kw)
        self._cache['coord_opts'] = _coord_opts

        return self._cache['coord']

    def get_error_samples(self, size=1, rnd=None):
        """Generate a sampling from the Gaia error distribution for each source.

        This function constructs the astrometric covariance matrix for each
        source and generates a specified number of random samples from the error
        distribution for each source. This does not handle spatially-dependent
        correlations. Samplings generated with this method can be used to, e.g.,
        propagate the Gaia errors through coordinate transformations or
        analyses.

        Parameters
        ----------
        size : int
            The number of random samples per soure to generate.
        rnd : ``numpy.random.RandomState``, optional
            The random state.

        Returns
        -------
        g_samples : `pyia.GaiaData`
            The same data table, but now each Gaia coordinate entry contains
            samples from the error distribution.

        """
        if rnd is None:
            rnd = np.random.RandomState()

        C = self.get_cov().copy()
        rv_mask = ~np.isfinite(C[:, 5, 5])
        C[rv_mask, 5, 5] = 0.

        arrs = []
        for k, unit in self._cache['cov_units'].items():
            arrs.append(getattr(self, k).to_value(unit))
        y = np.stack(arrs).T

        samples = np.array([rnd.multivariate_normal(y[i], C[i], size=size)
                            for i in range(len(y))])

        d = self.data.copy()
        for i, (k, unit) in enumerate(self._cache['cov_units'].items()):
            d[k] = samples[..., i] * unit

        return self.__class__(d)
