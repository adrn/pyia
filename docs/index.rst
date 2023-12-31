.. _pyia:

###################################################
pyia: a Python package for working with *Gaia* data
###################################################

Source code `on GitHub <https://github.com/adrn/pyia>`_.

:mod:`pyia` provides a class for working with data from the `ESA *Gaia* mission
<https://www.esa.int/Science_Exploration/Space_Science/Gaia>`_. Here is a short list of
features provided by this package:

* Access to *Gaia* data columns as Astropy :class:`~astropy.units.Quantity` objects with
  units associated (e.g., ``data.parallax`` will have associated units
  :attr:`~astropy.units.milliarcsecond`).
* Reconstruct covariance matrices for *Gaia* sources given errors and correlation
  coefficients (:meth:`pyia.GaiaData.get_cov`).
* Generate random samples from the astrometric error distribution for each source (a
  multivariate Gaussian centered on the catalog values)
  (:meth:`pyia.GaiaData.get_error_samples`).
* Create Astropy :class:`~astropy.coordinates.SkyCoord` objects from data
  (:attr:`pyia.GaiaData.skycoord`).
* Execute small remote queries via the *Gaia* science archive and automatically fetch
  the results (:meth:`pyia.GaiaData.from_query`).

************
Installation
************

Install :mod:`pyia` with ``pip``. We recommend installing the latest version from GitHub
directly with::

    pip install git+https://github.com/adrn/pyia

To install the latest stable version instead, use::

    pip install pyia


***************
Getting started
***************

The key class in this package is :class:`~pyia.GaiaData`. Creating an instance of this
class gives you access to the features listed above. :class:`~pyia.GaiaData` objects
can be created either by passing in a string filename, by passing in a pre-loaded
:class:`~astropy.table.Table` or :class:`~pandas.DataFrame` object, or by executing a
remote query with the :meth:`~pyia.GaiaData.from_query` class method. Let's now execute
some imports we will need later::

    >>> import astropy.table as at
    >>> import astropy.units as u
    >>> import numpy as np
    >>> import pyia


Creating a GaiaData instance
============================

From a filename
---------------

If you already have a file on disk with *Gaia* catalog data, you can create a new
:class:`~pyia.GaiaData` instance by passing in the path to the file::

    >>> g = pyia.GaiaData('path/to/file.fits')  # doctest: +SKIP

As an example, we have provided a small subset of *Gaia* DR3 data with this package in the ``tests/data`` directory. We can load this data with::

    >>> g = pyia.GaiaData('tests/data/gdr3_sm.fits')
    >>> g
    <GaiaData: 1000 rows>

The file can be in any format readable by Astropy, including FITS, HDF5, ASCII, CSV,
etc.


From a Table or DataFrame instance
----------------------------------

As mentioned above, you can also pass in a pre-loaded :class:`~astropy.table.Table` or
:class:`~pandas.DataFrame` object. For example, let's create a small table with some
data values::

    >>> tbl = at.Table(
    ...     {"ra": np.arange(10) * u.deg, "dec": np.linspace(-50, 50, 10) * u.deg}
    ... )
    >>> tbl  # doctest: +ELLIPSIS
    <Table length=10>
    ra           dec
    deg           deg
    ...

We could then pass this table instance directly in to :class:`~pyia.GaiaData`::

    >>> pyia.GaiaData(tbl)
    <GaiaData: 10 rows>


From an ADQL query to the Gaia archive
--------------------------------------

You can also create a :class:`~pyia.GaiaData` instance by executing a remote query to
the `*Gaia* archive <https://gea.esac.esa.int/archive/>`_ using the
:meth:`~pyia.GaiaData.from_query` class method. This method takes in a string ADQL query
and returns a new :class:`~pyia.GaiaData` instance with the results of the query. For
example, to query the *Gaia* archive for the 100 brightest sources within 10 degrees of
the Galactic anticenter::

    >>> query = """SELECT TOP 100 * FROM gaiadr3.gaia_source AS gaia
    ...     WHERE CONTAINS(
    ...         POINT(gaia.l, gaia.b),
    ...         CIRCLE(180, 0, 10)
    ...     )=1
    ...     AND phot_g_mean_mag IS NOT NULL
    ...     ORDER BY gaia.phot_g_mean_mag ASC
    ... """
    >>> pyia.GaiaData.from_query(query)  # doctest: +SKIP
    <GaiaData: 100 rows>


From a ``source_id`` from a different data release
--------------------------------------------------

Finally, you can create a :class:`~pyia.GaiaData` instance from a ``source_id`` from a
different *Gaia* data release. For example, to create a :class:`~pyia.GaiaData` instance
with data from *Gaia* DR3 given a ``source_id`` from *Gaia* DR2::

    >>> pyia.GaiaData.from_source_id(
    ...     4049398731303340416, source_id_dr="dr2", data_dr="dr3"
    ... )  # doctest: +SKIP
    <GaiaData: 1 rows>


Accessing data columns from a GaiaData instance and filtering values
====================================================================

For the subsequent examples below, we will use the small, random subset of *Gaia* DR3
data that is provided with this package::

    >>> g = pyia.GaiaData('tests/data/gdr3_sm.fits')

Once you have a :class:`~pyia.GaiaData` instance — named ``g`` here — you can access any
of the *Gaia* data columns as attributes. For example, to access the parallax and proper
motion in right ascension::

    >>> g.parallax[:4] # doctest: +FLOAT_CMP
    <Quantity [1.08492386, 0.85435642, 1.04200764, 0.58765983] mas>
    >>> g.pmra[:4] # doctest: +FLOAT_CMP
    <Quantity [ 2.54902781, -2.60680601, -1.72970386, -3.5660617 ] mas / yr>

The data is stored internally as an Astropy :class:`~astropy.table.Table` instance. You
can access this table directly using the ``.data`` attribute::

    >>> type(g.data)
    <class 'astropy.table.table.Table'>

The :class:`~pyia.GaiaData` object supports indexing and slicing like normal Python
array-like objects::

    >>> g[:4]
    <GaiaData: 4 rows>
    >>> g[g.parallax < 0.5*u.mas]
    <GaiaData: 653 rows>

You can also filter the data into rectangular selection regions over any set of columns
using the :meth:`~pyia.GaiaData.filter` method. This method takes keyword arguments
equal to column names, and an iterable specifying the minimum and maximum range to
select each data column within (pass ``None`` to allow any min or max). For example, to
select only sources with parallax less than 0.5 mas and proper motion in right ascension
between 5 mas/yr and 10 mas/yr::

    >>> g.filter(parallax=(None, 0.5*u.mas), pmra=[5, 10] * u.mas/u.yr)
    <GaiaData: 10 rows>


Dealing with negative parallaxes
================================

A large number of sources in *Gaia* have unmeasured parallaxes, some of which will be
negative. For example, in our sample data from *Gaia* DR3 loaded above, we can check how
many missing or negative parallaxes there are::

    >>> (~np.isfinite(g.parallax) | (g.parallax <= 0)).sum()
    245

For these sources, trying to transform the parallaxes to a distance will throw an error
because a negative parallax has no corresponding distance::

    >>> g.get_distance()  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    ValueError: some parallaxes are negative, which are not interpretable as distances. See the discussion in this paper: https://arxiv.org/abs/1507.02105 . You can convert negative parallaxes to NaN distances by providing the `allow_negative=True` argument.

We may therefore want to fill those values with NaN's or just filter them out when
computing the distances corresponding to the *Gaia* parallaxes. We can do this using the
:meth:`pyia.GaiaData.get_distance` method by passing in ``allow_negative=True``, which
by default sets the distances for any missing or negative parallaxes to NaN::

    >>> dist = g.get_distance(allow_negative=True)

For example, the first 8 values::

    >>> dist[:8]  # doctest: +FLOAT_CMP
    <Distance [ 921.72366887, 1170.47168706,  959.68586146, 1701.66471259,
                         nan, 2288.64173352,  420.24514217, 2839.42980486] pc>

We also provide a short-hand attribute for this method, ``.distance``, which also
replaces invalid parallax values with NaN values::

    >>> g.distance[:8]  # doctest: +FLOAT_CMP
    <Distance [ 921.72366887, 1170.47168706,  959.68586146, 1701.66471259,
                         nan, 2288.64173352,  420.24514217, 2839.42980486] pc>

We could alternatively set a real-valued fill value for these sources::

    >>> dist = g.get_distance(allow_negative=True, fill_value=1*u.megaparsec)
    >>> dist[:8]  # doctest: +FLOAT_CMP
    <Distance [9.21723669e+02, 1.17047169e+03, 9.59685861e+02, 1.70166471e+03,
               1.00000000e+06, 2.28864173e+03, 4.20245142e+02, 2.83942980e+03] pc>

Or specify a minimum parallax below which to set the distance to the fill value::

    >>> dist = g.get_distance(allow_negative=True, min_parallax=1*u.mas)
    >>> dist[:8]  # doctest: +FLOAT_CMP
    <Distance [921.72366887,          nan, 959.68586146,          nan,
                        nan,          nan, 420.24514217,          nan] pc>


Getting a SkyCoord object
=========================

We can retrieve a :class:`~astropy.coordinates.SkyCoord` object using the
:attr:`pyia.GaiaData.skycoord` attribute::

    >>> c = g.skycoord

Here, any missing or invalid distance (parallax) or radial velocity values will be
filled with NaN values. For example, the first 4 rows of our sample data::

    >>> c[:4] # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, pc)
        [(286.7169129 ,   0.2761946 ,  921.72366887),
         (152.62506771, -64.26769943, 1170.47168706),
         (348.43498221,  43.9407042 ,  959.68586146),
         (252.72971694, -37.91721471, 1701.66471259)]
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        [( 2.54902781, -4.07554383, nan), (-2.60680601,  4.41433169, nan),
         (-1.72970386, -3.35328852, nan), (-3.5660617 , -4.13250825, nan)]>

We could instead fill any invalid parallax and radial velocity measurements with real
values by passing in custom distance and/or radial velocity data to ``get_skycoord()``::

    >>> c = g.get_skycoord(
    ...     distance=g.get_distance(
    ...         min_parallax=1e-3 * u.mas,
    ...         fill_value=1 * u.megaparsec,
    ...         allow_negative=True
    ...     ),
    ...     radial_velocity=g.get_radial_velocity(fill_value=1e6 * u.km/u.s)
    ... )
    >>> c.galactic[:4]  # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b, distance) in (deg, deg, pc)
        [( 34.93906404,  -3.30633332,  921.72366887),
         (286.53235191,  -6.64869788, 1170.47168706),
         (104.93519089, -15.49052125,  959.68586146),
         (346.09665657,   4.14148012, 1701.66471259)]
     (pm_l_cosb, pm_b, radial_velocity) in (mas / yr, mas / yr, km / s)
        [(-2.46002534, -4.12987598, 1000000.),
         (-4.68065797,  2.09122042, 1000000.),
         (-2.87915693, -2.43862148, 1000000.),
         (-5.45701325,  0.12420493, 1000000.)]>


Generating error samples
========================

*Gaia* provides estimates of the uncertainties on the astrometric measurements for each
coordinate along with correlation coefficients between the different measurements. These
are stored in the, e.g., ``.parallax_error``, ``.pmra_error``, ``.pmdec_error``, etc.
and, e.g., ``.ra_parallax_corr``, ``.parallax_pmra_corr``, etc. columns. We can use
these to reconstruct the full covariance matrix for each source, and then generate
random samples from the error distribution for each source. This can be useful when, for
example, transforming to a new coordinate system when you want to fully propagate the
uncertainties in the *Gaia* data through your analysis.

To retrieve the covariance matrix itself, use the :meth:`pyia.GaiaData.get_cov` method::

    >>> cov, cov_units = g.get_cov()
    >>> cov.shape
    (1000, 6, 6)

This returns a 3D array with shape ``(N, 6, 6)``, where ``N`` is the number of sources.
By default, the returned units are: (deg, deg, milliarcsecond, milliarcsecond/year,
milliarcsecond/year, km/s), and these units are also returned by the
:meth:`pyia.GaiaData.get_cov` method::

    >>> cov_units
    {'ra': Unit("deg"), 'dec': Unit("deg"), 'parallax': Unit("mas"), 'pmra': Unit("mas / yr"), 'pmdec': Unit("mas / yr"), 'radial_velocity': Unit("km / s")}

The units can be changed in the resulting covariance matrix by passing in a dictionary
with the new desired units. For example, to get the RA and Dec terms in radians instead
of degrees::

    >>> cov2 = g.get_cov(units=dict(ra=u.radian, dec=u.radian))

To generate random samples from the *Gaia* error distribution for each source, use the
:meth:`pyia.GaiaData.get_error_samples` method. As a demonstration, we will generate
32 error samples for sources in our subset data set that have well-measured parallaxes
and non-null radial velocities::

    >>> g_subset = g.filter(
    ...     parallax_over_error=(8, None), radial_velocity=[-1000, 1000] * u.km/u.s
    ... )
    >>> g_subset_samples = g_subset.get_error_samples(size=32)

To see the full set of features provided by the :class:`~pyia.GaiaData` class, see the
API documentation below.


***
API
***

.. automodapi:: pyia
    :no-inheritance-diagram:
