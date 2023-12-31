# Standard library
import pathlib
import warnings

import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import pytest

# Third-party
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename

# Project
from pyia.data import GaiaData, gaia_unit_map

dr_filenames = [
    get_pkg_data_filename("data/gdr2_sm.fits"),
    get_pkg_data_filename("data/gedr3_sm.fits"),
    get_pkg_data_filename("data/gdr3_sm.fits"),
]


@pytest.mark.parametrize(
    ("filename", "fmt"),
    [
        (get_pkg_data_filename("data/gdr2_sm.fits"), "fits"),
        (get_pkg_data_filename("data/gdr2_sm.vot"), "votable"),
    ],
)
def test_init(filename, fmt):
    """Test initializing the GaiaData class with different options"""

    if any("fits" in x for x in pathlib.Path(filename).suffixes):
        gd = GaiaData(fits.getdata(filename, 1))
        assert isinstance(gd.data, Table)

        kw = {"unit_parse_strict": "silent"}
    else:
        kw = {}

    gd = GaiaData(Table.read(filename, format=fmt, **kw))
    gd = GaiaData(filename, format=fmt)


@pytest.mark.parametrize("filename", dr_filenames)
def test_slicing_getattr(filename):
    gd = GaiaData(filename)

    for k, unit in gaia_unit_map.items():
        if k in gd.data.columns:
            # have to special case the flux columns:
            if unit == u.ph / u.s:
                continue
            assert getattr(gd, k).unit == unit

    size = len(gd)
    assert size > 0

    gd_slc = gd[:4]
    assert len(gd_slc) == 4

    gd_one = gd[4]
    assert len(gd_one) == 1

    assert "parallax" in dir(gd)
    assert "source_id" in dir(gd)
    gd.source_id  # noqa: B018
    gd.parallax  # noqa: B018
    gd.pmra  # noqa: B018

    gd["source_id"]
    gd["parallax"]


@pytest.mark.parametrize("filename", dr_filenames)
def test_str_repr(filename):
    gd = GaiaData(filename)
    assert "GaiaData" in repr(gd)
    assert "1000 rows" in repr(gd)


@pytest.mark.parametrize("filename", dr_filenames)
def test_computed_quantities(filename):
    gd = GaiaData(filename)

    pm = gd.get_pm()
    assert pm.shape == (len(gd), 2)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert isinstance(gd.distance, coord.Distance)
        assert gd.distance.shape == (len(gd),)

        assert gd.distmod.unit == u.mag

        assert gd.vtan.unit == u.km / u.s
        assert gd.vtan.shape == (len(gd), 2)


@pytest.mark.parametrize("filename", dr_filenames)
def test_cov(filename):
    kw = {}
    if "gedr3" in filename:
        kw["radial_velocity_colname"] = "dr2_radial_velocity"
        kw["radial_velocity_error_colname"] = "dr2_radial_velocity_error"
    gd = GaiaData(filename, **kw)

    C, _ = gd.get_cov()
    assert C.shape == (len(gd), 6, 6)

    with pytest.raises(RuntimeError):
        gd.get_cov(RAM_threshold=1 * u.kilobyte)  # force failure


@pytest.mark.parametrize("filename", dr_filenames)
def test_skycoord(filename):
    kw = {}
    if "gedr3" in filename:
        kw["radial_velocity_colname"] = "dr2_radial_velocity"
        kw["radial_velocity_error_colname"] = "dr2_radial_velocity_error"
    gd = GaiaData(filename, **kw)
    gd = gd[np.isfinite(gd.get_radial_velocity())]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = gd.skycoord
        assert len(c) == len(gd)

        c = gd.get_skycoord()
        assert np.allclose(
            gd.get_radial_velocity(),
            c.radial_velocity,
            equal_nan=True,
            atol=1e-12,
            rtol=0,
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gd["dist"] = (
            coord.Distance(parallax=gd.parallax, allow_negative=True).kpc * u.kpc
        )
        c = gd.get_skycoord(distance="dist")
        assert np.allclose(
            gd.get_radial_velocity(),
            c.radial_velocity,
            equal_nan=True,
            atol=1e-12,
            rtol=0,
        )
        assert np.allclose(gd.distance, c.distance, equal_nan=True, atol=1e-12, rtol=0)


@pytest.mark.parametrize("filename", dr_filenames)
def test_setattr(filename):
    tbl = Table.read(filename, unit_parse_strict="silent")
    tbl["arr_column"] = np.arange(len(tbl)) * 10.0
    gd = GaiaData(tbl)

    # Setting a quantity column
    with pytest.raises(ValueError, match="Quantity-like object"):
        gd.parallax = np.arange(len(gd))

    new_vals = np.arange(len(gd)) * u.microarcsecond
    gd.parallax = new_vals
    assert np.all(np.asarray(gd.data["parallax"]) == new_vals.value)
    assert gd.parallax.unit == u.microarcsecond

    # Setting an array column
    new_vals = np.arange(len(gd))
    gd.arr_column = new_vals
    assert np.all(gd.data["arr_column"] == new_vals)


@pytest.mark.remote_data()
def test_from_query():
    q = """SELECT TOP 10 * FROM gaiadr1.tgas_source"""
    gd = GaiaData.from_query(q)

    assert len(gd) == 10


@pytest.mark.remote_data()
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
@pytest.mark.parametrize(
    ("filename", "source_dr"),
    [
        ("data/gdr2_sm.fits", "dr2"),
        ("data/gdr3_sm.fits", "dr3"),
        ("data/gedr3_sm.fits", "edr3"),
    ],
)
def test_from_source_id(filename, source_dr):
    filename = get_pkg_data_filename(filename)
    source_id = Table.read(filename)["source_id"][0]

    for data_dr in ["dr2", "edr3", "dr3"]:
        if (source_dr == "dr3" and data_dr == "edr3") or (
            source_dr == "edr3" and data_dr == "dr3"
        ):
            continue

        gd = GaiaData.from_source_id(source_id, source_id_dr=source_dr, data_dr=data_dr)
        assert len(gd) > 0
        assert gd.designation[0].startswith(f"Gaia {data_dr.upper()}")


@pytest.mark.parametrize("filename", dr_filenames)
def test_get_samples(filename):
    kw = {}
    if "gedr3" in filename:
        kw["radial_velocity_colname"] = "dr2_radial_velocity"
        kw["radial_velocity_error_colname"] = "dr2_radial_velocity_error"
    gd = GaiaData(filename, **kw)

    g_samples = gd.get_error_samples(size=16)
    assert g_samples.ra.shape == (len(gd), 16)
    assert g_samples.dec.shape == (len(gd), 16)
    assert g_samples.parallax.shape == (len(gd), 16)
    assert g_samples.pmra.shape == (len(gd), 16)
    assert g_samples.pmdec.shape == (len(gd), 16)
    assert g_samples.get_radial_velocity().shape == (len(gd), 16)

    c = g_samples.get_skycoord(distance=False)
    assert c.shape == (1000, 16)


@pytest.mark.parametrize("filename", dr_filenames)
def test_extinction(filename):
    g = GaiaData(filename)

    rng = np.random.default_rng(seed=42)
    ebv = rng.uniform(0, 0.1, len(g))
    g.get_ext(ebv=ebv)
    g.get_BP0()
    g.get_BP0(ebv=ebv)


@pytest.mark.parametrize("filename", dr_filenames)
def test_filter(filename):
    g = GaiaData(filename)

    new_g = g.filter(parallax=(0.5, 5) * u.mas, phot_g_mean_mag=(None, 17.5 * u.mag))

    assert new_g.parallax.min() > 0.5 * u.mas
    assert new_g.parallax.max() < 5 * u.mas
    assert new_g.phot_g_mean_mag.max() < 17.5 * u.mag


class TestCustomRVDist:
    def setup_class(self):
        # First, make some renamed columns:
        filename = get_pkg_data_filename("data/gdr3_sm.fits")
        tbl = Table.read(filename, unit_parse_strict="silent")
        self.orig_tbl = tbl.copy()

        tbl = tbl[~tbl["radial_velocity"].mask & (tbl["parallax"] > 0)]

        tbl["VHELIO"] = tbl["radial_velocity"]
        tbl["VERR"] = tbl["radial_velocity_error"]

        rng = np.random.default_rng(seed=42)
        plx_samples = rng.normal(
            tbl["parallax"], tbl["parallax_error"], size=(128, len(tbl))
        ).T
        dd = 1 / plx_samples
        tbl["dist50"] = np.nanmedian(dd, axis=1) * u.kpc
        tbl["dist_err"] = (
            np.diff(np.nanpercentile(dd, [16, 84], axis=1), axis=0)[0] / 2 * u.kpc
        )
        tbl.remove_columns(
            ["parallax", "parallax_error", "radial_velocity", "radial_velocity_error"]
        )
        self.tbl = tbl

        self.kw = {
            "radial_velocity_colname": "VHELIO",
            "radial_velocity_error_colname": "VERR",
            "distance_colname": "dist50",
            "distance_error_colname": "dist_err",
        }

    def test_success_init(self):
        g = GaiaData(self.tbl, **self.kw)

        d = g.distance
        assert np.all(d.value == self.tbl["dist50"])

        rv = g.get_radial_velocity()
        assert np.all(rv.value == self.tbl["VHELIO"])

        # Pass unit in explicitly:
        tbl = self.tbl.copy()
        tbl["VHELIO"].unit = None
        tbl["dist50"].unit = None
        g = GaiaData(tbl, **self.kw, radial_velocity_unit=u.m / u.s, distance_unit=u.pc)
        assert g.get_radial_velocity().unit == u.m / u.s
        assert g.get_distance().unit == u.pc

    def test_fail_init(self):
        # Colname that is missing:
        for k in self.kw:
            kk = self.kw.copy()
            kk[k] = "not_a_column"
            with pytest.raises(ValueError, match="not found in data table"):
                GaiaData(self.tbl, **kk)

        # No associated unit:
        for v in self.kw.values():
            tbl = self.tbl.copy()
            tbl[v].unit = None
            with pytest.raises(ValueError, match="does not have a unit"):
                GaiaData(tbl, **self.kw)

    def test_skycoord(self):
        g = GaiaData(self.tbl, **self.kw)
        c = g.get_skycoord()
        assert u.allclose(c.radial_velocity, g.VHELIO)
        assert u.allclose(c.distance, g.dist50)

        c = g.get_skycoord(radial_velocity=g.VHELIO * 2)
        assert not u.allclose(c.radial_velocity, g.VHELIO)

    def test_cov(self):
        g = GaiaData(self.tbl, **self.kw)
        C, units = g.get_cov()
        assert C.shape == (len(g), 6, 6)
        assert np.allclose(C[:, 2, 2], g.dist_err.to_value(g.dist50.unit) ** 2)
        assert units["dist50"] == g.dist50.unit
        for i in range(6):
            if i == 2:
                continue
            assert np.allclose(C[:, 2, i], 0.0)
            assert np.allclose(C[:, i, 2], 0.0)

        # Default behavior warns if missing correlation coefficients
        with pytest.warns(RuntimeWarning, match="Missing correlation"):
            g.get_cov(warn_missing_corr=True)

        # Subset of coordinates
        C, units = g.get_cov(coords=["ra", "dec", "dist50"], warn_missing_corr=False)
        assert C.shape == (len(g), 3, 3)

    def test_error_samples(self):
        g = GaiaData(self.tbl, **self.kw)
        s = g.get_error_samples(size=16, warn_missing_corr=False)
        assert s.ra.shape == (len(g), 16)
        assert s.dec.shape == (len(g), 16)
        assert s.VHELIO.shape == (len(g), 16)
        assert s.dist50.shape == (len(g), 16)

        g = GaiaData(self.orig_tbl)
        s = g.get_error_samples(size=16, coords=["ra", "dec", "parallax"])
        assert s.ra.shape == (len(g), 16)
        assert s.parallax.shape == (len(g), 16)
        assert s.pmra.shape == (len(g),)
