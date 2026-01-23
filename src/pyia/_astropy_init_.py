# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import warnings

from astropy.config.configuration import (
    ConfigurationDefaultMissingError,
    ConfigurationDefaultMissingWarning,
)

try:
    from astropy.config.configuration import update_default_config
except ImportError:  # Astropy >= 6
    update_default_config = None

config_dir = os.path.dirname(__file__)

if update_default_config:
    try:
        update_default_config(__package__, config_dir)
    except ConfigurationDefaultMissingError as e:
        wmsg = (
            e.args[0]
            + " Cannot install default profile. If you are importing from source, this is expected."
        )
        warnings.warn(wmsg, ConfigurationDefaultMissingWarning)