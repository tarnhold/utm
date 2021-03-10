from setuptools import setup, find_packages
from Cython.Build import cythonize

from utm._version import __version__


setup(
    name='utm',
    version=__version__,
    author='Tobias Bieniek',
    author_email='Tobias.Bieniek@gmx.de',
    url='https://github.com/Turbo87/utm',
    description='Bidirectional UTM-WGS84 converter for python',
    keywords=['utm', 'wgs84', 'coordinate', 'converter'],
    classifiers=[
        'Programming Language :: Python',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Development Status :: 4 - Beta',
        'Environment :: Other Environment',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: GIS',
    ],
    packages=['utm'],
    ext_modules=cythonize('utm/conversion.pyx', compiler_directives={'language_level': 3}),
    scripts=['scripts/utm-converter'],
    test_suite='test',
)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
