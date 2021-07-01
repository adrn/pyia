1.3 (unreleased)
================

- Added a new ``GaiaData.from_source_id()`` classmethod to enable retrieving
  Gaia data from some data release, given a source_id from another data release.

- Allow ``GaiaData`` instances to be created from a ``pathlib.Path`` object


1.2 (2020-12-01)
================

- Added support for EDR3.


1.1 (2020-07-13)
================

- Fix packaging error.

1.0 (2020-07-11)
================

- Added ``get_uwe`` and ``get_ruwe`` methods to compute the unit weight error
  and renormalized unit weight error (using code from Anthony Brown).

0.4 (2019-03-24)
================

- Added ``get_distance`` and ``get_radial_velocity`` methods for filling
  unmeasured or invalid values.
- Added support for computing extinction values and extinction-corrected
  photometry for the Gaia DR2 bandpasses.
- Added a method for generating random samples from the Gaia astrometric error
  distribution for each source (``get_error_samples()``).
- Output ``SkyCoord`` objects now contain an ``obstime`` corresponding to the
  reference epoch of the data.

0.3 (2018-06-04)
================

- Changed underlying data structure to be an Astropy Table instead of a Pandas
  DataFrame.
- Fixed handling columns that have units or unrecognized units.

0.2 (2018-04-24)
================

- More efficient loading of FITS files
- Added all units for all expected columns

0.1 (2018-04-22)
================

- First release
