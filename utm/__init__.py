from utm.conversion import to_latlon, from_latlon, latlon_to_zone_number, latitude_to_zone_letter, mod_angle
from utm.error import OutOfRangeError

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "0.0.0+unknown"

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
