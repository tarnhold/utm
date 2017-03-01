import math
from utm.error import OutOfRangeError

__all__ = ['to_latlon', 'from_latlon']

K0 = 0.9996

# GRS80
R = 6378137.0
F = 1.0 / 298.257222101

# WGS84
#R = 6378137.0
#F = 1.0 / 298.257223563

# first numerical eccentricity squared (e^2)
E = 2.0 * F - math.pow(F, 2.0) # [-]
E2 = E ** 2
E3 = E ** 3
E4 = E ** 4
E5 = E ** 5
# second numerical eccentricity squared (e^'2)
E_P2 = E / (1.0 - E)

SQRT_E = math.sqrt(1 - E)
_E = (1 - SQRT_E) / (1 + SQRT_E)
_E2 = _E ** 2
_E3 = _E ** 3
_E4 = _E ** 4
_E5 = _E ** 5

# Meridian distance from latitude
#
# Deakin, R. E. (2006): Meridian Distance, School of Mathematical &
#   Geospatial Sciences, RMIT University, Melbourne, March 2006.
# Weintritt, A. (2013): So, What is Actually the Distance from the
#   Equator to the Pole? - Overview of Meridian Distance Approximations,
#   TransNav the International Journal on Marine Navigation and Safety
#   of Sea Transportation, Volume 7, Number 2, June 2013.
M1 =                  (1 - E / 4 - 3 * E2 / 64 -  5 * E3 / 256 - 175 * E4 / 16384 - 441 * E5 / 65536)
M2 =   3. /       8 * (    E +         E2 /  4 + 15 * E3 / 128 +  35 * E4 /   512 + 735 * E5 / 16384)
M3 =  15. /     256 * (                E2      +  3 * E3 /   4 +  35 * E4 /    64 + 105 * E5 /   256)
M4 =  35. /    3072 * (                               E3       +   5 * E4 /     4 + 315 * E5 /   256)
M5 = 315. /  131072 * (                                                E4         +   7 * E5 /     4)
M6 = 693. / 1310720 *                                                                     E5

# Latitude from distance
#
# Deakin, R. E. (2012): Great Elliptic Arc Distance, School of Mathematical &
#   Geospatial Sciences, RMIT University, Melbourne, January 2012.
P2 = (   3. /    2 * _E  -  27. /  32 * _E3 + 269. / 512 * _E5)
P3 = (  21. /   16 * _E2 -  55. /  32 * _E4)
P4 = ( 151. /   96 * _E3 - 417. / 128 * _E5)
P5 = (1097. /  512 * _E4)
P6 = (8011. / 2560 * _E5) # 13th place difference

R = 6378137

ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"


def to_latlon(easting, northing, zone_number, zone_letter=None, northern=None, strict=True):
    """This function convert an UTM coordinate into Latitude and Longitude

        Parameters
        ----------
        easting: int
            Easting value of UTM coordinate

        northing: int
            Northing value of UTM coordinate

        zone number: int
            Zone Number is represented with global map numbers of an UTM Zone
            Numbers Map. More information see utmzones [1]_

        zone_letter: str
            Zone Letter can be represented as string values. Where UTM Zone
            Designators can be accessed in [1]_

        northern: bool
            You can set True or False to set this parameter. Default is None


       .. _[1]: http://www.jaworski.ca/utmzones.htm

    """
    if not zone_letter and northern is None:
        raise ValueError('either zone_letter or northern needs to be set')

    elif zone_letter and northern is not None:
        raise ValueError('set either zone_letter or northern, but not both')

    if strict:
        if not 100000 <= easting < 1000000:
            raise OutOfRangeError('easting out of range (must be between 100.000 m and 999.999 m)')
        if not 0 <= northing <= 10000000:
            raise OutOfRangeError('northing out of range (must be between 0 m and 10.000.000 m)')
    if not 1 <= zone_number <= 60:
        raise OutOfRangeError('zone number out of range (must be between 1 and 60)')

    if zone_letter:
        zone_letter = zone_letter.upper()

        if not 'C' <= zone_letter <= 'X' or zone_letter in ['I', 'O']:
            raise OutOfRangeError('zone letter out of range (must be between C and X)')

        northern = (zone_letter >= 'N')

    x = easting - 500000
    y = northing

    if not northern:
        y -= 10000000

    m = y / K0
    mu = m / (R * M1)

    # meridian distance to latitude
    p_rad = (mu +
             P2 * math.sin(2 * mu) +
             P3 * math.sin(4 * mu) +
             P4 * math.sin(6 * mu) +
             P5 * math.sin(8 * mu) +
             P6 * math.sin(10 * mu))

    p_sin2 = math.sin(p_rad) ** 2
    p_cos = math.cos(p_rad)

    p_tan = math.tan(p_rad)
    p_tan2 = p_tan ** 2
    p_tan4 = p_tan ** 4
    p_tan6 = p_tan ** 6

    n = R / math.sqrt(1 - E * p_sin2)
    r = (1 - E) / (1 - E * p_sin2)

    c = _E * p_cos**2
    c2 = c ** 2
    c3 = c ** 3
    c4 = c ** 4

    d = x / (n * K0)

    # Kelly, Kevin M. (1986): Coordinate Transformations - Universal
    #   Transverse Mercator/Geographic. Ontario Ministry of Natural
    #   Resources. March 1986.
    # Snyder, John P. (1987): Map Projections - A Working Manual, U.S.
    #   Geological Survey Professional Paper 1395, p.60ff, Washington.
    # Hofmann-Wellenhof, B.; Kienast, G.; Lichtenegger, H. (1994): GPS in der
    #   Praxis, Springer-Verlag Wien New York, p.97ff, Wien.
    latitude = (p_rad - (p_tan / r) *
                (d ** 2 / 2 -
                 d ** 4 / 24 * (5 + 3 * p_tan2 + c - 4 * c2 - 9 * p_tan2 * c) +
                 d ** 6 / 720 * (61 + 90 * p_tan2 + 45 * p_tan4 + 46 * c - 3 * c2 + 100 * c3 + 88 * c4 - 252 * p_tan2 * c - 66 * p_tan2 * c2 + 84 * p_tan4 * c3 - 192 * p_tan2 * c4 - 90 * p_tan4 * c + 225 * p_tan4 * c2) -
                 d ** 8 / 40320 * (1385 + 3633 * p_tan2 + 4095 * p_tan4 + 1575 * p_tan6)))

    longitude = (d -
                 d ** 3 / 6 * (1 + 2 * p_tan2 + c) +
                 d ** 5 / 120 * (5 + 28 * p_tan2 + 24 * p_tan4 + 6 * c - 3 * c2 - 4 * c3 + 8 * p_tan2 * c + 4 * p_tan2 * c2 + 24 * p_tan2 * c3) -
                 d ** 7 / 5040 * (61 + 662 * p_tan2 + 1320 * p_tan4 + 720 * p_tan6)) / p_cos

    return (math.degrees(latitude),
            math.degrees(longitude) + zone_number_to_central_longitude(zone_number))


def from_latlon(latitude, longitude, force_zone_number=None):
    """This function convert Latitude and Longitude to UTM coordinate

        Parameters
        ----------
        latitude: float
            Latitude between 80 deg S and 84 deg N, e.g. (-80.0 to 84.0)

        longitude: float
            Longitude between 180 deg W and 180 deg E, e.g. (-180.0 to 180.0).

        force_zone number: int
            Zone Number is represented with global map numbers of an UTM Zone
            Numbers Map. You may force conversion including one UTM Zone Number.
            More information see utmzones [1]_

       .. _[1]: http://www.jaworski.ca/utmzones.htm
    """
    if not -80.0 <= latitude <= 84.0:
        raise OutOfRangeError('latitude out of range (must be between 80 deg S and 84 deg N)')
    if not -180.0 <= longitude <= 180.0:
        raise OutOfRangeError('longitude out of range (must be between 180 deg W and 180 deg E)')

    lat_rad = math.radians(latitude)

    lat_tan = math.tan(lat_rad)
    lat_tan2 = lat_tan ** 2
    lat_tan4 = lat_tan ** 4
    lat_tan6 = lat_tan ** 6

    if force_zone_number is None:
        zone_number = latlon_to_zone_number(latitude, longitude)
    else:
        zone_number = force_zone_number

    zone_letter = latitude_to_zone_letter(latitude)

    lon_rad = math.radians(longitude)
    central_lon = zone_number_to_central_longitude(zone_number)
    central_lon_rad = math.radians(central_lon)

    n = R / math.sqrt(1 - E * sin(lat_rad)**2)
    c = E_P2 * cos(lat_rad)**2
    c2 = c ** 2
    c3 = c ** 3
    c4 = c ** 4

    a = cos(lat_rad) * (lon_rad - central_lon_rad)

    # meridian distance from latitude
    m = R * (M1 * lat_rad -
             M2 * math.sin(2 * lat_rad) +
             M3 * math.sin(4 * lat_rad) -
             M4 * math.sin(6 * lat_rad) +
             M5 * math.sin(8 * lat_rad) -
             M6 * math.sin(10 * lat_rad))

    # Kelly, Kevin M. (1986): Coordinate Transformations - Universal
    #    Transverse Mercator/Geographic. Ontario Ministry of Natural
    #    Resources. March 1986.
    # Snyder, John P. (1987): Map Projections - A Working Manual, U.S.
    #   Geological Survey Professional Paper 1395, p.60ff, Washington.
    # Hofmann-Wellenhof, B.; Kienast, G.; Lichtenegger, H. (1994): GPS in der
    #    Praxis, Springer-Verlag Wien New York, p.97ff, Wien.
    easting = K0 * n * (a +
                        a ** 3 / 6 * (1 - lat_tan2 + c) +
                        a ** 5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 14 * c + 13 * c2 + 4 * c3 - 58 * c * lat_tan2 - 64 * c2 * lat_tan2 - 24 * lat_tan2 * c3 ) +
                        a ** 7 / 5040 * (61 - 479 * lat_tan2 + 179 * lat_tan4 - lat_tan6)
                        ) + 500000

    northing = K0 * (m + n * lat_tan * (a ** 2 / 2 +
                                        a ** 4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c2) +
                                        a ** 6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 270 * c  + 445 * c2 + 324 * c3 + 88 * c4 - 330 * c * lat_tan2 - 680 * lat_tan2 * c2 - 600 * lat_tan2 * c3 - 192 * lat_tan2 * c4) +
                                        a ** 8 / 40320 * (1385 - 3111 * lat_tan2 + 543 * lat_tan4 - lat_tan6)
                                        ))

    if latitude < 0:
        northing += 10000000

    return easting, northing, zone_number, zone_letter


def latitude_to_zone_letter(latitude):
    if -80 <= latitude <= 84:
        return ZONE_LETTERS[int(latitude + 80) >> 3]
    else:
        return None


def latlon_to_zone_number(latitude, longitude):
    if 56 <= latitude < 64 and 3 <= longitude < 12:
        return 32

    if 72 <= latitude <= 84 and longitude >= 0:
        if longitude <= 9:
            return 31
        elif longitude <= 21:
            return 33
        elif longitude <= 33:
            return 35
        elif longitude <= 42:
            return 37

    return int((longitude + 180) / 6) + 1


def zone_number_to_central_longitude(zone_number):
    return (zone_number - 1) * 6 - 180 + 3

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
