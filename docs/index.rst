.. _pyia:

###################################################
pyia: a Python package for working with *Gaia* data
###################################################

What ``pyia`` can do for you:

* Provides access to Gaia data columns as `~astropy.units.Quantity` objects,
  i.e. with units (e.g., ``data.parallax`` will have units
  'milliarcsecond'),
* Construct covariance matrices for Gaia data (`pyia.GaiaData.get_cov()`),
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
example, I've created a small subset of `Gaia DR2 mock data
<http://dc.zah.uni-heidelberg.de/__system__/dc_tables/show/tableinfo/gdr2mock.main>`_
to use for examples below::

    >>> import astropy.units as u
    >>> from pyia import GaiaData
    >>> g = GaiaData('docs/_static/gdr2mock.fits')
    >>> g
    <GaiaData: 100 rows>

With this object, we can access any of the Gaia table column names using
attribute access. By using the attribute, we get back an
`~astropy.units.Quantity` object with the correct units::

    >>> g.parallax[:4]
    <Quantity [0.167184 , 0.134295 , 0.23638  , 0.0976651] marcsec>
    >>> g.phot_g_mean_mag[:4]
    <Quantity [20.6775, 20.4755, 19.6778, 20.4051] mag>

To access the raw data (stored internally as a Pandas `~pandas.DataFrame`), you can use the ``.data`` attribute::

    >>> type(g.data)
    <class 'pandas.core.frame.DataFrame'>

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
    >>> cov
    array([[[2.20150978e-13,            nan,            nan,            nan,
                    nan, 0.00000000e+00],
        [           nan, 1.44516248e-13,            nan,            nan,
                    nan, 0.00000000e+00],
        [           nan,            nan, 4.61578083e+00,            nan,
                    nan, 0.00000000e+00],
        [           nan,            nan,            nan, 1.43958020e+00,
                    nan, 0.00000000e+00],
        [           nan,            nan,            nan,            nan,
         9.68711734e-01, 0.00000000e+00],
        [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
         0.00000000e+00,            inf]],

       [[1.54148517e-13,            nan,            nan,            nan,
                    nan, 0.00000000e+00],
        [           nan, 1.01189510e-13,            nan,            nan,
                    nan, 0.00000000e+00],
        [           nan,            nan, 3.23194528e+00,            nan,
                    nan, 0.00000000e+00],
        [           nan,            nan,            nan, 1.00798643e+00,
                    nan, 0.00000000e+00],
        [           nan,            nan,            nan,            nan,
         6.78286731e-01, 0.00000000e+00],
        [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
         0.00000000e+00,            inf]]])

By default, the returned units are: (deg, deg, milliarcsecond,
milliarcsecond/year, milliarcsecond/year, km/s), but these can be changed by
passing in a dictionary with the new desired units::

    >>> cov2 = g[:2].get_cov(units=dict(ra=u.radian, dec=u.radian))
    >>> cov[0, 0, 0] # doctest: +FLOAT_CMP
    2.201509780681496e-13
    >>> cov2[0, 0, 0] # doctest: +FLOAT_CMP
    6.70618285130544e-17

We can also retrieve other useful quantities from this object, like a 2D proper
motion array with units (as a `~astropy.units.Quantity`)::

    >>> g[:4].pm
    <Quantity [[-0.26951  ,  0.0829121],
               [-0.618642 , -5.09182  ],
               [-3.29806  , -1.31502  ],
               [-4.84104  , -8.83496  ]] marcsec / yr>

Finally, we can retrieve a `~astropy.coordinates.SkyCoord` object for all rows::

    >>> c = g.skycoord
    >>> c[:4]
    <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, pc)
        [(265.74811032, -36.35797425,  5981.43359375),
         (265.74947102, -36.36186489,  7446.29394531),
         (265.74951857, -36.36185265,  4230.4765625 ),
         (265.75036833, -36.36044967, 10239.07226562)]
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        [(-0.26951 ,  0.0829121, -134.197  ),
         (-0.618642, -5.09182  , -233.08   ),
         (-3.29806 , -1.31502  ,    4.73733),
         (-4.84104 , -8.83496  ,  -63.567  )]>

But note that this computes the distance using 1/parallax.

***
API
***

.. automodapi:: pyia
    :no-inheritance-diagram:
    :skip: test
