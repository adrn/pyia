# Standard library
import abc
from collections import namedtuple

__all__ = ['ComponentMapping', 'DataMeta']


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
