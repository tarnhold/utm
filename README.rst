utm
===

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://github.com/Turbo87/utm/blob/master/LICENSE

Geodetic UTM/GRS80 projection

Modified version of the original UTM/WGS84 projection from https://github.com/Turbo87/utm

Because there was only a rough UTM projection formula in use, there was a
need to further expand the power series for geodetic precision (like this
is realized with pyproj or GeographicLib). Formulas are cross-checked and
based on geodetic literature. 

Modifications

* Expanded power series from meter to sub-millimeter precise projection
* Projection onto GRS80 ellipsoid, as this is the basis for ETRS89 coordinates
* Cythonization for performance
* Extended testing, invoked by setup.py test
* Some minor bugfixes
