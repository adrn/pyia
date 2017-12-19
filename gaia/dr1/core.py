# Standard library
import abc
from collections import namedtuple

# Third-party
from astropy.utils.compat.misc import override__dir__
import astropy.coordinates as coord
import astropy.units as u
import numpy as np

__all__ = ['ComponentMapping', 'DataMeta', 'BaseData', 'BaseGaiaData']


_ComponentMappingBase = \
    namedtuple('ComponentMapping',
               ('obj_name', 'tbl_name', 'tbl_unit'))


class ComponentMapping(_ComponentMappingBase):
    """
    This `~collections.namedtuple` is used to map from raw column names to
    common attribute names. For example, from ``pmra`` to ``pm_ra_cosdec``.
    """

    def __new__(cls, obj_name, tbl_name, tbl_unit):
        # this trick just provides some defaults
        return super().__new__(cls, obj_name, tbl_name, tbl_unit)


class DataMeta(abc.ABCMeta):
    def __new__(mcls, name, bases, members):

        if 'component_mappings' in members:
            mappings = members['component_mappings']
            found_mappings = True

        else:
            found_mappings = False

        # somewhat hacky, but this is the best way to get the MRO according to
        # https://mail.python.org/pipermail/python-list/2002-December/167861.html
        tmp_cls = super().__new__(mcls, name, bases, members)

        for m in (c.__dict__ for c in tmp_cls.__mro__):
            if (not found_mappings and 'component_mappings' in m):
                mappings = m['component_mappings']
                found_mappings = True

            if found_mappings:
                break
        else:
            raise ValueError("Couldn't find component_mappings frame "
                             "attribute.")

        members['_tbl_to_obj'] = dict([(m.tbl_name, m.obj_name)
                                       for m in mappings])
        members['_obj_to_tbl'] = dict([(m.obj_name, m.tbl_name)
                                       for m in mappings])
        members['_tbl_to_unit'] = dict([(m.tbl_name, m.tbl_unit)
                                        for m in mappings])
        members['_obj_to_unit'] = dict([(m.obj_name, m.tbl_unit)
                                        for m in mappings])

        # Cache a list of extra names to add to dir() so autocomplete works on
        # the object attribute names instead of the table column names
        members['_extra_dir_names'] = list(members['_obj_to_tbl'].keys())

        # Pre-compute unit scale factors to make generating the covariance
        # matrix fast. This is for cases where the '<name>_error' unit is
        # different from the unit of the column '<name>'
        _err_unit_scale_factor = dict()
        for name in members['_tbl_to_unit']:
            err_name = "{0}_error".format(name)
            if (name in members['_tbl_to_unit'] and
                    err_name in members['_tbl_to_unit']):
                err_unit = members['_tbl_to_unit'][err_name]
                col_unit = members['_tbl_to_unit'][name]
                _err_unit_scale_factor[err_name] = err_unit.to(col_unit)
        members['_err_unit_scale_factor'] = _err_unit_scale_factor

        return super().__new__(mcls, name, bases, members)


class BaseData(metaclass=DataMeta):

    component_mappings = []

    # --------------------------------------------------------------------------
    # Python special methods:
    #
    def __getattr__(self, name):
        # to prevent recursion errors:
        # http://nedbatchelder.com/blog/201010/surprising_getattr_recursion.html
        if name == 'data':
            raise AttributeError()

        if name in self._obj_to_tbl:
            tbl_name = self._obj_to_tbl[name]
            unit = self._obj_to_unit[name]
            return self.data[tbl_name] * unit

        elif name in self.data.colnames:
            return self.data[name]

        else:
            raise AttributeError("Object {} has no attribute '{}' and source "
                                 "data table has no column with that name."
                                 .format(self, name))

    def __len__(self):
        return len(self.data)

    @override__dir__
    def __dir__(self):
        """Override the builtin `dir` behavior to include representation names.
        """
        colnames = []
        for name in self.data.colnames:
            if (name in self._tbl_to_obj and
                    self._tbl_to_obj[name] not in self._extra_dir_names):
                colnames.append(name)

            elif name not in self._tbl_to_obj:
                colnames.append(name)

        return self._extra_dir_names + colnames


# Both TGAS and the Gaia source tables should work with the below
class BaseGaiaData(BaseData):

    # --------------------------------------------------------------------------
    # Other convenience methods:
    #
    def get_distance(self, lutz_kelker=False):
        """Return the distance with or without the Lutz-Kelker correction.

        Parameters
        ----------
        lutz_kelker : bool, optional
            Apply the Lutz-Kelker correction to the distance. You probably don't
            want to do this.
        """

        if lutz_kelker:
            snr = self.parallax / self.parallax_error
            tmp = self.parallax * (0.5 + 0.5*np.sqrt(1 - 16/snr**2))

        else:
            tmp = self.parallax

        return coord.Distance(tmp.to(u.pc, u.parallax()))

    def get_vtan(self, lutz_kelker=False):
        """Return the tangential velocity computed using the proper motion and
        distance.

        Parameters
        ----------
        lutz_kelker : bool, optional
            Apply the Lutz-Kelker correction to the distance. You probably don't
            want to do this.
        """
        d = self.get_distance(lutz_kelker=lutz_kelker)
        vra = (self.pm_ra_cosdec * d).to(u.km/u.s,
                                         u.dimensionless_angles()).value
        vdec = (self.pm_dec * d).to(u.km/u.s, u.dimensionless_angles()).value
        return np.vstack((vra, vdec)) * u.km/u.s

    def get_coord(self, lutz_kelker=False):
        """
        Return an `astropy.coordinates` object to represent all coordinates.

        Parameters
        ----------
        lutz_kelker : bool, optional
            Apply the Lutz-Kelker correction to the distance. You probably don't
            want to do this.
        """
        return coord.ICRS(ra=self.ra, dec=self.dec,
                          distance=self.get_distance(lutz_kelker=lutz_kelker),
                          pm_ra_cosdec=self.pm_ra_cosdec,
                          pm_dec=self.pm_dec)

    # --------------------------------------------------------------------------
    # Custom attributes:
    #
    @property
    def parallax_snr(self):
        """Parallax signal-to-noise ratio."""
        return self.parallax / self.parallax_error

    @property
    def distance(self):
        """Distance, computed as 1/parallax."""
        return self.get_distance(lutz_kelker=False)
