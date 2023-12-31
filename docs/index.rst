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

TODO: distance, skycoord, absolute mag?

Finally, we can retrieve a `~astropy.coordinates.SkyCoord` object for all rows::

    >>> c = g.skycoord
    >>> c[:4] # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, pc)
        [(192.37472723, -58.59207336, 4352.66598258),
         (218.43888361, -76.53529493, 1409.831102  ),
         (328.7882822 ,  52.06426886, 2110.71574001),
         (104.61548808,   2.60168493, 1188.72069818)]
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        [( -7.08493408,  1.49887161, nan), (-10.43896891, -7.55166654, nan),
         ( -2.67527741, -1.83343139, nan), ( -2.84067421,  4.59017005, nan)]>

But note that this computes the distance using 1/parallax.


Dealing with negative parallaxes or custom distances
====================================================

A large number of sources in *Gaia* have unmeasured parallaxes, some of which will
even be negative. For these sources, naively transforming to an astropy
coordinate object will throw an error (because a negative parallax has no
corresponding distance!). We may want to fill those values with NaN's or just
filter them out. Let's work now with a small subset of the *Gaia* data that
contains some negative parallax measurements::

    >>> g = pyia.GaiaData(data_path / 'gdr2_sm_negplx.fits')
    >>> len(g)
    8
    >>> g[np.isfinite(g.parallax) & (g.parallax > 0)]
    <GaiaData: 4 rows>

So of the 8 sources, 4 have parallax values that physically can't be converted
into distance values. If we try to get a coordinate object for the whole table,
we will therefore get an error from ``astropy.coordinates``::

    >>> g.get_skycoord() # doctest: +SKIP
    ...
    ERROR: ValueError: Some parallaxes are negative, which are notinterpretable as distances. See the discussion in this paper: https://arxiv.org/abs/1507.02105 . If you want parallaxes to pass through, with negative parallaxes instead becoming NaN, use the `allow_negative=True` argument. [astropy.coordinates.distances]

To get distance values for this table and fill the unmeasured parallaxes, we can
instead use the ``.get_distance()`` method::

    >>> dist = g.get_distance(min_parallax=1e-3*u.mas)
    >>> dist # doctest: +FLOAT_CMP
    <Distance [          nan,           nan, 2408.32836149,           nan,
               4310.30352187,           nan, 1273.15031694,  873.80934793] pc>

We can then pass in our filled distance values to ``get_skycoord()`` to retrieve
a coordinate object with any invalid distance values filled::

    >>> c = g.get_skycoord(distance=dist)

If we fill with NaN values, this will protect us from accidentally interpreting
the distance values as real with subsequent analysis. However, the NaN distances
will cause any coordinate transformations to fail (the transformation will fill
all values with NaN when a distance is NaN)::

    >>> c = g.get_skycoord(distance=dist)
    >>> c.galactic # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b, distance) in (deg, deg, pc)
        [(         nan,        nan,           nan),
         (         nan,        nan,           nan),
         (288.32933722, 2.46100214, 2408.32836149),
         (         nan,        nan,           nan),
         (288.2192576 , 2.45677656, 4310.30352187),
         (         nan,        nan,           nan),
         (288.32950975, 2.50348803, 1273.15031694),
         (288.03472562, 1.89930315,  873.80934793)]
     (pm_l_cosb, pm_b, radial_velocity) in (mas / yr, mas / yr, km / s)
        [(nan, nan, nan), (nan, nan, nan), (nan, nan, nan), (nan, nan, nan),
         (nan, nan, nan), (nan, nan, nan), (nan, nan, nan), (nan, nan, nan)]>

We can therefore fill any invalid parallax and radial velocity measurements by
passing in custom distance and/or radial velocity data to ``get_skycoord()``::

    >>> c = g.get_skycoord(distance=g.get_distance(min_parallax=1e-3*u.mas,
    ...                                            parallax_fill_value=1e-5*u.mas),
    ...                    radial_velocity=g.get_radial_velocity(fill_value=1e8*u.km/u.s))
    >>> c.galactic # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b, distance) in (deg, deg, pc)
        [(288.19236138, 2.43659588, 1.00000000e+08),
         (288.32526751, 2.18391781, 1.00000000e+08),
         (288.32933722, 2.46100214, 2.40832836e+03),
         (288.24762666, 2.24721944, 1.00000000e+08),
         (288.2192576 , 2.45677656, 4.31030352e+03),
         (288.38223571, 2.50727897, 1.00000000e+08),
         (288.32950975, 2.50348803, 1.27315032e+03),
         (288.03472562, 1.89930315, 8.73809348e+02)]
     (pm_l_cosb, pm_b, radial_velocity) in (mas / yr, mas / yr, km / s)
        [(         nan,         nan,    nan),
         ( -6.49472182,  0.75047358, 1.e+08),
         ( -3.87451798, -0.8264891 , 1.e+08),
         ( -6.52193605, -3.52263043, 1.e+08),
         ( -3.34309802, -0.92741245, 1.e+08),
         ( -5.49737417, -4.03517439, 1.e+08),
         ( -7.35557765, -1.27131657, 1.e+08),
         (-12.82935235, -1.37765851, 1.e+08)]>


Generating error samples
========================

TODO: cov matrix and error samples

We can also get the 6 by 6 covariance matrices for all rows in the current data
object. As a demo, we'll just get the covariance matrix for the first 2 rows
(note that in the mock DR2 sample, all off-diagonal terms are set to nan, so
only the diagonal elements exist!)::

    >>> cov = g[:2].get_cov()
    >>> cov.shape
    (2, 6, 6)
    >>> cov # doctest: +SKIP
    array([[[ 1.69220176e-17,  3.21641380e-18, -2.95728307e-11,
             -5.06610201e-11, -1.42767067e-11,  0.00000000e+00],
            [ 3.21641380e-18,  2.04162196e-17, -3.07238505e-11,
             -9.54818920e-12, -6.23362213e-11,  0.00000000e+00],
            [-2.95728307e-11, -3.07238505e-11,  5.27824144e-04,
              1.78662319e-04,  2.42514756e-04,  0.00000000e+00],
            [-5.06610201e-11, -9.54818920e-12,  1.78662319e-04,
              7.26667008e-04,  2.80791812e-04,  0.00000000e+00],
            [-1.42767067e-11, -6.23362213e-11,  2.42514756e-04,
              2.80791812e-04,  7.50643258e-04,  0.00000000e+00],
            [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              0.00000000e+00,  0.00000000e+00,             inf]],

           [[ 2.88454811e-16,  2.23388900e-17, -1.86116307e-10,
             -7.52130802e-10, -6.85621463e-11,  0.00000000e+00],
            [ 2.23388900e-17,  2.19942577e-16, -1.86746169e-10,
              9.53999432e-11, -3.54438368e-10,  0.00000000e+00],
            [-1.86116307e-10, -1.86746169e-10,  5.03113212e-03,
             -6.74448280e-04,  2.56024316e-03,  0.00000000e+00],
            [-7.52130802e-10,  9.53999432e-11, -6.74448280e-04,
              1.28494693e-02,  1.76042822e-03,  0.00000000e+00],
            [-6.85621463e-11, -3.54438368e-10,  2.56024316e-03,
              1.76042822e-03,  1.15835427e-02,  0.00000000e+00],
            [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              0.00000000e+00,  0.00000000e+00,             inf]]])

By default, the returned units are: (deg, deg, milliarcsecond,
milliarcsecond/year, milliarcsecond/year, km/s), but these can be changed by
passing in a dictionary with the new desired units::

    >>> cov2 = g[:2].get_cov(units=dict(ra=u.radian, dec=u.radian))
    >>> cov[0, 0, 2] # doctest: +FLOAT_CMP
    -2.957283073704542e-11
    >>> cov2[0, 0, 2] # doctest: +FLOAT_CMP
    -5.16143265496424e-13

We can also retrieve other useful quantities from this object, like a 2D proper
motion array with units (as a `~astropy.units.Quantity`)::

    >>> g[:4].pm # doctest: +FLOAT_CMP
    <Quantity [[ -7.08493408,   1.49887161],
               [-10.43896891,  -7.55166654],
               [ -2.67527741,  -1.83343139],
               [ -2.84067421,   4.59017005]] marcsec / yr>




It is sometimes useful to generate random samples from the *Gaia* error
distribution for each source. This can be useful when, for example, transforming
to a new coordinate system when you want to propagate the (correlated!)
uncertainty in the *Gaia* data through your analysis. We can generate samples from
the *Gaia* error distribution using ``pyia``. As an example, we'll work with a
small subset of the *Gaia* data that have radial velocity measurements, sub-select
only nearby sources, and then generate error samples for the sources. We'll then
transform the samples to Galactocentric coordinates to look at the uncertainty
distribution for the full-space velocity.

First, let's load the data::

    >>> g_rv = pyia.GaiaData(data_path / 'gdr2_rv_sm.fits')

All of these sources have measured radial velocities::

    >>> g_rv.radial_velocity[:4] # doctest: +FLOAT_CMP
    <Quantity [  7.89796709,  30.88496542,   3.04709697, -34.91701273] km / s>

Let's now select only nearby (within 500 pc) sources::

    >>> g_rv = g_rv[g_rv.parallax > 2*u.mas]

To generate samples from the error distribution, we use the
``.get_error_samples()`` method, and pass in the number of samples to generate
(here, 256)::

    >>> import numpy as np
    >>> g_samples = g_rv.get_error_samples(size=256,
    ...                                    rnd=np.random.RandomState(seed=42))

Let's now get a ``SkyCoord`` object to represent the data for these sources and
samples, and transform to a Galactocentric coordinate frame using the Astropy
coordinate transformation machinery::

    >>> c_samples = g_samples.get_skycoord()
    >>> import astropy.coordinates as coord
    >>> _ = coord.galactocentric_frame_defaults.set('v4.0')
    >>> galcen = c_samples.transform_to(coord.Galactocentric)

Let's now look at the uncertainty on the magnitude of the total velocity, ``v``,
for each of these sources::

    >>> v = galcen.velocity.norm()

And finally, let's compute (from the error samples) the uncertainty on the total
velocity for these sources::

    >>> err_v = np.std(v, axis=1)
    >>> err_v # doctest: +FLOAT_CMP
    <Quantity [1.19963698, 0.90448355, 0.76143172, 0.58234578, 0.56352556,
               5.05012474, 0.57285428, 0.41450777, 0.13843949, 0.40697535,
               1.01139842, 1.72859784, 0.10139791, 2.64303227, 0.22751807,
               0.0897801 , 0.62391201, 0.65292275, 0.77665343, 0.84242905,
               1.27174119, 2.13547719, 1.80547001, 1.15228657] km / s>

Most of these uncertainties are less than 1-2 km/s! These take into account the
parallax, proper motion, and radial velocity uncertainties provided by Gaia.


***
API
***

.. automodapi:: pyia
    :no-inheritance-diagram:
