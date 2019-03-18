.. _pyia:

###################################################
pyia: a Python package for working with *Gaia* data
###################################################

Source code `on GitHub <https://github.com/adrn/pyia>`_.

What ``pyia`` can do for you:

* Provides access to Gaia data columns as `~astropy.units.Quantity` objects,
  i.e. with units (e.g., ``data.parallax`` will have units
  'milliarcsecond'),
* Construct covariance matrices for Gaia data (`pyia.GaiaData.get_cov()`),
* Generate random samples from the Gaia error distribution per source
  (`pyia.GaiaData.get_error_samples()`),
* Create `~astropy.coordinates.SkyCoord` objects from Gaia data
  (`pyia.GaiaData.skycoord`),
* Support for executing simple (small) remote queries via the Gaia science
  archive and automatically fetching results (`pyia.GaiaData.from_query()`).

************
Installation
************

See instructions here:

.. toctree::
   :maxdepth: 1

   install

***************
Getting started
***************

The key class in this package is the `~pyia.GaiaData` class. Creating an
instance of this class gives you access to the features mentioned above.
`~pyia.GaiaData` objects can be created either by passing in a string filename,
by passing in a pre-loaded `~astropy.table.Table` or `~pandas.DataFrame` object,
or by executing a remote query with `~pyia.GaiaData.from_query()`.

If you've already downloaded some Gaia data in tabular format and saved it to a
file on disk, you can create an instance by passing the path to the file. As an
example, I've created a small subset of Gaia DR2 data downloaded from the
Gaia science archive to use for examples below::

    >>> import astropy.units as u
    >>> from pyia import GaiaData
    >>> g = GaiaData('docs/_static/gdr2_sm.fits')
    >>> g
    <GaiaData: 100 rows>

As mentioned above, you can also pass in pre-loaded data::

    >>> from astropy.table import Table
    >>> tbl = Table.read('docs/_static/gdr2_sm.fits')
    >>> GaiaData(tbl)
    <GaiaData: 100 rows>

With this object, we can access any of the Gaia table column names using
attribute access. By using the attribute, we get back an
`~astropy.units.Quantity` object with the correct units::

    >>> g.parallax[:4] # doctest: +FLOAT_CMP
    <Quantity [0.37942248, 0.12635312, 0.29292803, 0.62822535] mas>
    >>> g.phot_g_mean_mag[:4] # doctest: +FLOAT_CMP
    <Quantity [18.38548 , 19.645763, 18.432325, 19.573648] mag>

To access the raw data (stored internally as an Astropy `~astropy.table.Table`), you can use the ``.data`` attribute::

    >>> type(g.data)
    <class 'astropy.table.table.Table'>

The `~pyia.GaiaData` object supports indexing and slicing like normal Python lists / Numpy arrays / Astropy tables::

    >>> g[:4]
    <GaiaData: 4 rows>
    >>> g[g.parallax < 0.1*u.mas]
    <GaiaData: 7 rows>

We can also get the 6 by 6 covariance matrices for all rows in the current data
object. As a demo, we'll just get the covariance matrix for the first 2 rows
(note that in the mock DR2 sample, all off-diagonal terms are set to nan, so
only the diagonal elements exist!)::

    >>> cov = g[:2].get_cov()
    >>> cov.shape
    (2, 6, 6)
    >>> cov # doctest: +SKIP
    array([[[ 1.11300727e-15,  3.39476922e-17,  6.32569989e-10,
             -1.83063606e-09, -3.14384117e-10,  0.00000000e+00],
            [ 3.39476922e-17,  1.24216618e-15,  1.47796850e-09,
             -3.26525998e-10, -3.24502195e-09,  0.00000000e+00],
            [ 6.32569989e-10,  1.47796850e-09,  2.06969264e-02,
             -3.80907519e-03, -3.41596832e-03,  0.00000000e+00],
            [-1.83063606e-09, -3.26525998e-10, -3.80907519e-03,
              6.94569085e-02, -9.45256223e-03,  0.00000000e+00],
            [-3.14384117e-10, -3.24502195e-09, -3.41596832e-03,
             -9.45256223e-03,  5.71660391e-02,  0.00000000e+00],
            [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              0.00000000e+00,  0.00000000e+00,             inf]],

           [[ 7.09961079e-15,  5.03257689e-16,  5.48053482e-09,
             -1.17483063e-08, -2.40896552e-09,  0.00000000e+00],
            [ 5.03257689e-16,  1.06594795e-14,  1.24732958e-08,
             -3.91322912e-11, -3.50557367e-08,  0.00000000e+00],
            [ 5.48053482e-09,  1.24732958e-08,  1.60209695e-01,
              2.21121108e-02,  2.88825197e-02,  0.00000000e+00],
            [-1.17483063e-08, -3.91322912e-11,  2.21121108e-02,
              5.30656499e-01, -6.01228583e-02,  0.00000000e+00],
            [-2.40896552e-09, -3.50557367e-08,  2.88825197e-02,
             -6.01228583e-02,  5.71332641e-01,  0.00000000e+00],
            [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              0.00000000e+00,  0.00000000e+00,             inf]]])

By default, the returned units are: (deg, deg, milliarcsecond,
milliarcsecond/year, milliarcsecond/year, km/s), but these can be changed by
passing in a dictionary with the new desired units::

    >>> cov2 = g[:2].get_cov(units=dict(ra=u.radian, dec=u.radian))
    >>> cov[0, 0, 2] # doctest: +FLOAT_CMP
    2.201509780681496e-13
    >>> cov2[0, 0, 2] # doctest: +FLOAT_CMP
    6.70618285130544e-17

We can also retrieve other useful quantities from this object, like a 2D proper
motion array with units (as a `~astropy.units.Quantity`)::

    >>> g[:4].pm # doctest: +FLOAT_CMP
    <Quantity [[-4.24880717, -4.52180435],
               [ 1.90317677, -3.80754094],
               [-2.89919833, -3.87763849],
               [-4.62431593, -2.00130269]] mas / yr>

Finally, we can retrieve a `~astropy.coordinates.SkyCoord` object for all rows::

    >>> c = g.skycoord
    >>> c[:4] # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, pc)
        [(321.33207659, 50.18344643, 2635.58449057),
         (322.35324537, 50.63947336, 7914.32760316),
         (322.13767285, 50.36279294, 3413.80786286),
         (322.10360875, 50.62762822, 1591.7854934 )]
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        [(-4.24880717, -4.52180435, nan), ( 1.90317677, -3.80754094, nan),
         (-2.89919833, -3.87763849, nan), (-4.62431593, -2.00130269, nan)]>

But note that this computes the distance using 1/parallax.


Generating error samples
========================

It is sometimes useful to generate random samples from the Gaia error
distribution for each source. This can be useful when, for example, transforming
to a new coordinate system when you want to propagate the (correlated!)
uncertainty in the Gaia data through your analysis. We can generate samples from
the Gaia error distribution using ``pyia``. As an example, we'll work with a
small subset of the Gaia data that have radial velocity measurements, sub-select
only nearby sources, and then generate error samples for the sources. We'll then
transform the samples to Galactocentric coordinates to look at the uncertainty
distribution for the full-space velocity.

First, let's load the data:

    >>> g_rv = GaiaData('docs/_static/gdr2_rv_sm.fits')

All of these sources have measured radial velocities:

    >>> g_rv.radial_velocity[:4] # doctest: +FLOAT_CMP
    <Quantity [  7.89796709,  30.88496542,   3.04709697, -34.91701273] km / s>

Let's now select only nearby (within 500 pc) sources:

    >>> g_rv = g_rv[g_rv.parallax > 2*u.mas]

To generate samples from the error distribution, we use the
``.get_error_samples()`` method, and pass in the number of samples to generate
(here, 256):

    >>> g_samples = g_rv.get_error_samples(size=256)

Let's now get a ``SkyCoord`` object to represent the data for these sources and
samples, and transform to a Galactocentric coordinate frame using the Astropy
coordinate transformation machinery:

    >>> c_samples = g_samples.get_skycoord()
    >>> import astropy.coordinates as coord
    >>> galcen = c_samples.transform_to(coord.Galactocentric)

Let's now look at the uncertainty on the magnitude of the total velocity, ``v``,
for each of these sources:

    >>> v = galcen.velocity.norm()

And finally, let's compute (from the error samples) the uncertainty on the total
velocity for these sources:


    >>> err_v = np.std(v, axis=1)
    >>> err_v # doctest: +FLOAT_CMP
    <Quantity [1.34468499, 0.84924219, 0.8439712 , 0.51272844, 0.56964407,
               5.09773664, 0.49509594, 0.38900574, 0.12777623, 0.33003514,
               0.90980165, 1.55980192, 0.10561152, 2.24834014, 0.23978129,
               0.08849996, 0.60940975, 0.66246121, 0.72935829, 0.73371712,
               1.39479966, 1.98871529, 1.67979077, 1.22117557] km / s>

Most of these uncertainties are less than 1-2 km/s! These take into account the
parallax, proper motion, and radial velocity uncertainties provided by Gaia.


***
API
***

.. automodapi:: pyia
    :no-inheritance-diagram:
    :skip: test
