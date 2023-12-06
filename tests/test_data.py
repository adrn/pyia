# Standard library
from __future__ import annotations

import pathlib

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
    get_pkg_data_filename("data/gdr3_sm.fits"),
]


@pytest.mark.parametrize(
    "filename,fmt",
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
    gd.source_id
    gd.parallax
    gd.pmra

    gd["source_id"]
    gd["parallax"]


@pytest.mark.parametrize("filename", dr_filenames)
def test_str_repr(filename):
    gd = GaiaData(filename)
    assert "GaiaData" in repr(gd)
    assert "100 rows" in repr(gd)


@pytest.mark.parametrize("filename", dr_filenames)
def test_computed_quantities(filename):
    gd = GaiaData(filename)

    pm = gd.get_pm()
    assert pm.shape == (len(gd), 2)

    assert isinstance(gd.distance, coord.Distance)
    assert gd.distance.shape == (len(gd),)

    assert gd.distmod.unit == u.mag

    assert gd.vtan.unit == u.km / u.s
    assert gd.vtan.shape == (len(gd), 2)


@pytest.mark.parametrize("filename", dr_filenames)
def test_cov(filename):
    gd = GaiaData(filename)

    C, _ = gd.get_cov()
    assert C.shape == (len(gd), 6, 6)

    with pytest.raises(RuntimeError):
        gd.get_cov(RAM_threshold=1 * u.kilobyte)  # force failure


@pytest.mark.parametrize("filename", dr_filenames)
def test_skycoord(filename):
    gd = GaiaData(filename)
    gd = gd[np.isfinite(gd.radial_velocity)]

    c = gd.skycoord
    assert len(c) == len(gd)

    c = gd.get_skycoord(radial_velocity="radial_velocity")
    assert np.all(c.radial_velocity == gd.radial_velocity)

    gd["dist"] = coord.Distance(parallax=gd.parallax).kpc * u.kpc
    c = gd.get_skycoord(radial_velocity="radial_velocity", distance="dist")
    assert np.all(c.radial_velocity == gd.radial_velocity)
    assert np.all(c.distance == gd.distance)


@pytest.mark.parametrize("filename", dr_filenames)
def test_setattr(filename):
    tbl = filename
    tbl["arr_column"] = np.arange(len(tbl)) * 10.0
    gd = GaiaData(tbl)

    # Setting a quantity column
    with pytest.raises(ValueError):
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
def test_from_source_id():
    filename = get_pkg_data_filename("data/gdr2_sm.fits")
    tbl = filename

    gd = GaiaData.from_source_id(tbl["source_id"][0], "dr2", "dr2")
    assert len(gd) == 1
    assert gd.designation[0].startswith("Gaia DR2")

    gd = GaiaData.from_source_id(tbl["source_id"][0], "dr2", data_dr="edr3")
    assert len(gd) == 1
    assert gd.designation[0].startswith("Gaia EDR3")


@pytest.mark.parametrize("filename", dr_filenames)
def test_get_samples(filename):
    gd = GaiaData(filename)

    g_samples = gd.get_error_samples(size=16)
    assert g_samples.ra.shape == (len(gd), 16)
    assert g_samples.dec.shape == (len(gd), 16)
    assert g_samples.parallax.shape == (len(gd), 16)
    assert g_samples.pmra.shape == (len(gd), 16)
    assert g_samples.pmdec.shape == (len(gd), 16)
    assert g_samples.radial_velocity.shape == (len(gd), 16)

    c = g_samples.get_skycoord(distance=False)
    assert c.shape == (100, 16)


def test_compute_ruwe():
    filename = get_pkg_data_filename("data/gdr2_sm.fits")
    gd = GaiaData(filename)
    ruwe = gd.get_ruwe()
    assert len(ruwe) == len(gd)


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
