import utm as UTM
import unittest


class UTMTestCase(unittest.TestCase):
    def assert_utm_equal(self, a, b, precision=6):
        self.assertAlmostEqual(a[0], b[0], precision)
        self.assertAlmostEqual(a[1], b[1], precision)
        self.assertEqual(a[2], b[2])
        self.assertEqual(a[3].upper(), b[3].upper())

    def assert_latlon_equal(self, a, b, precision=5):
        self.assertAlmostEqual(a[0], b[0], precision)
        self.assertAlmostEqual(a[1], b[1], precision)

class KnownValuesGRS80(UTMTestCase):
    # Known UTM values were projected from latitude and longitude values
    # using GeographicLib (onto GRS80 ellipsoid!). As this library has a
    # much higher series expansion and a different implementation we can
    # assume they are more accurate and use this as reference.
    known_values = [
        # Aachen, Germany
        (
            (50.77534556, 6.08388667),
            (294408.662941387, 5628897.512984829, 32, 'U'),
            {'northern': True},
        ),
        # New York, USA
        (
            (40.71435000, -74.00597000),
            (583959.959045332, 4507523.086854665, 18, 'T'),
            {'northern': True},
        ),
        # Wellington, New Zealand
        (
            (-41.28646000, 174.77623611),
            (313783.980049117, 5427057.313755062, 60, 'G'),
            {'northern': False},
        ),
        # Capetown, South Africa
        (
            (-33.92486889, 18.42405500),
            (261877.350976653, 6243185.700844696, 34, 'H'),
            {'northern': False},
        ),
        # Mendoza, Argentina
        (
            (-32.89018000, -68.84405000),
            (514586.227836383, 6360876.825073616, 19, 'h'),
            {'northern': False},
        ),
        # Fairbanks, Alaska, USA
        (
            (64.83777806, -147.71638889),
            (466013.322449279, 7190567.781669118, 6, 'W'),
            {'northern': True},
        ),
        # Ben Nevis, Scotland, UK
        (
            (56.79680000, -5.00601000),
            (377485.765670114, 6296561.854117111, 30, 'V'),
            {'northern': True},
        ),
        # Latitude 84
        (
            (84, -5.00601),
            (476594.34011230164, 9328501.361833721, 30, 'X'),
            {'northern': True},
        ),
    ]

    def test_from_latlon(self):
        '''from_latlon should give known result with known input'''
        for latlon, utm, _ in self.known_values:
            result = UTM.from_latlon(*latlon)
            self.assert_utm_equal(utm, result)

    def test_to_latlon(self):
        '''to_latlon should give known result with known input'''
        for latlon, utm, utm_kw in self.known_values:
            result = UTM.to_latlon(*utm)
            self.assert_latlon_equal(latlon, result)

            result = UTM.to_latlon(*utm[0:3], **utm_kw)
            self.assert_latlon_equal(latlon, result)

    def test_from_latlon_roundtrip(self):
        '''from_latlon look how good roundtrip fits'''
        for latlon, utm, utm_kw in self.known_values:
            utmr = UTM.from_latlon(*latlon)
            result = UTM.to_latlon(*utmr[0:3], **utm_kw)
            # we should get the same values as the initial input
            self.assert_latlon_equal(latlon, result, 5)

    # series expansion for inverse (utm->llh) is worse than forward (llh->utm)
    # calculation so we cannot expect same accuracy
    def test_to_latlon_roundtrip(self):
        '''to_latlon look how good roundtrip fits'''
        for latlon, utm, utm_kw in self.known_values:
            latlonr = UTM.to_latlon(*utm)

            # disable strict lat/lon range check, because roundtrip
            # of "Latitude 84" is 84.00000016... which is outside range
            result = UTM.from_latlon(*latlonr, strict=False)
            # we should get the same values as the initial input
            self.assert_latlon_equal(utm, result, 0)


class BadInput(UTMTestCase):
    def test_from_latlon_range_checks(self):
        '''from_latlon should fail with out-of-bounds input'''
        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, -100, 0)
        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, -80.1, 0)

        # test valid range
        for i in range(-8000, 8400):
            UTM.from_latlon(i / 100.0, 0)

        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, 84.1, 0)
        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, 100, 0)

        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, 0, -300)
        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, 0, -180.1)

        # test valid range
        for i in range(-18000, 18000):
            UTM.from_latlon(0, i / 100.0)

        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, 0, 180.1)
        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, 0, 300)

        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, -100, -300)
        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, 100, -300)
        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, -100, 300)
        self.assertRaises(UTM.OutOfRangeError, UTM.from_latlon, 100, 300)

    def test_to_latlon_range_checks(self):
        '''to_latlon should fail with out-of-bounds input'''

        # validate input
        self.assertRaises(
            ValueError, UTM.to_latlon, 500000, -100000, 32, 'U', northern=True)
        self.assertRaises(
            ValueError, UTM.to_latlon, 500000, -100000, 32, '')
        self.assert_latlon_equal((0.904730614584, 9.0),
            UTM.to_latlon(500000, 100000, 32, '', northern=True))

        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, -100000, 32, 'UU')


        # test easting range
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 0, 5000000, 32, 'U')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 99999, 5000000, 32, 'U')

        # valid range
        for i in range(100000, 999999, 1000):
            UTM.to_latlon(i, 5000000, 32, 'U')

        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 1000000, 5000000, 32, 'U')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 100000000000, 5000000, 32, 'U')

        # test northing range
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, -100000, 32, 'U')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, -1, 32, 'U')

        # valid range
        for i in range(10, 10000000, 1000):
            UTM.to_latlon(500000, i, 32, 'U')

        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 10000001, 32, 'U')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 50000000, 32, 'U')

        # test zone numbers
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, -1, 'U')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 0, 'U')

        # valid range
        for i in range(1, 60):
            UTM.to_latlon(500000, 5000000, i, 'U')

        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 61, 'U')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 1000, 'U')

        # test zone letters
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 32, 'A')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 32, 'B')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 32, 'I')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 32, 'O')

        # there are no zone numbers 32, 34 and 36 in X
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 32, 'X')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 34, 'X')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 36, 'X')

        # valid range
        for i in range(ord('C'), ord('X')):
            i = chr(i)
            if i != 'I' and i != 'O':
                UTM.to_latlon(500000, 5000000, 32, i)

        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 32, 'Y')
        self.assertRaises(
            UTM.OutOfRangeError, UTM.to_latlon, 500000, 5000000, 32, 'Z')


class SpecialZones(unittest.TestCase):

    def assert_zone_equal(self, result, expected_number, expected_letter):
        self.assertEqual(result[2], expected_number)
        self.assertEqual(result[3].upper(), expected_letter.upper())

    # test 31X, 33X, 35X, 37X
    def test_zones_X(self):
        # test lower left and upper left
        self.assert_zone_equal(UTM.from_latlon(72, 0), 31, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 9), 33, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 21), 35, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 33), 37, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 42), 38, 'X')

        self.assert_zone_equal(UTM.from_latlon(84, 0), 31, 'X')
        self.assert_zone_equal(UTM.from_latlon(84, 9), 33, 'X')
        self.assert_zone_equal(UTM.from_latlon(84, 21), 35, 'X')
        self.assert_zone_equal(UTM.from_latlon(84, 33), 37, 'X')
        self.assert_zone_equal(UTM.from_latlon(84, 42), 38, 'X')

        # test inside
        self.assert_zone_equal(UTM.from_latlon(72, 6), 31, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 12), 33, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 18), 33, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 24), 35, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 30), 35, 'X')
        self.assert_zone_equal(UTM.from_latlon(72, 36), 37, 'X')

    # test 31V and 32V
    def test_inside(self):
        # test 31V
        self.assert_zone_equal(UTM.from_latlon(56, 0), 31, 'V')
        self.assert_zone_equal(UTM.from_latlon(56, 2.999999), 31, 'V')

        # test 32V
        self.assert_zone_equal(UTM.from_latlon(56, 3), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(56, 6), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(56, 9), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(56, 11.999999), 32, 'V')

        self.assert_zone_equal(UTM.from_latlon(60, 3), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(60, 6), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(60, 9), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(60, 11.999999), 32, 'V')

        self.assert_zone_equal(UTM.from_latlon(63.999999, 3), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(63.999999, 6), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(63.999999, 9), 32, 'V')
        self.assert_zone_equal(UTM.from_latlon(63.999999, 11.999999), 32, 'V')

    def test_left_of(self):
        self.assert_zone_equal(UTM.from_latlon(55.999999, 2.999999), 31, 'U')
        self.assert_zone_equal(UTM.from_latlon(56, 2.999999), 31, 'V')
        self.assert_zone_equal(UTM.from_latlon(60, 2.999999), 31, 'V')
        self.assert_zone_equal(UTM.from_latlon(63.999999, 2.999999), 31, 'V')
        self.assert_zone_equal(UTM.from_latlon(64, 2.999999), 31, 'W')

    def test_right_of(self):
        self.assert_zone_equal(UTM.from_latlon(55.999999, 12), 33, 'U')
        self.assert_zone_equal(UTM.from_latlon(56, 12), 33, 'V')
        self.assert_zone_equal(UTM.from_latlon(60, 12), 33, 'V')
        self.assert_zone_equal(UTM.from_latlon(63.999999, 12), 33, 'V')
        self.assert_zone_equal(UTM.from_latlon(64, 12), 33, 'W')

    def test_below(self):
        self.assert_zone_equal(UTM.from_latlon(55.999999, 3), 31, 'U')
        self.assert_zone_equal(UTM.from_latlon(55.999999, 6), 32, 'U')
        self.assert_zone_equal(UTM.from_latlon(55.999999, 9), 32, 'U')
        self.assert_zone_equal(UTM.from_latlon(55.999999, 11.999999), 32, 'U')
        self.assert_zone_equal(UTM.from_latlon(55.999999, 12), 33, 'U')

    def test_above(self):
        self.assert_zone_equal(UTM.from_latlon(64, 3), 31, 'W')
        self.assert_zone_equal(UTM.from_latlon(64, 6), 32, 'W')
        self.assert_zone_equal(UTM.from_latlon(64, 9), 32, 'W')
        self.assert_zone_equal(UTM.from_latlon(64, 11.999999), 32, 'W')
        self.assert_zone_equal(UTM.from_latlon(64, 12), 33, 'W')


if __name__ == '__main__':
    unittest.main()

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
