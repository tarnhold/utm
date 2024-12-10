from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name='utm',
    use_scm_version={"write_to": "utm/_version.py"},
    author='Thomas Arnhold',
    author_email='thomas@arnhold.org',
    url='https://github.com/tarnhold/utm',
    description='Bidirectional UTM/GRS80 converter for python',
    keywords=['utm', 'grs80', 'coordinate', 'converter'],
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
    setup_requires=[
       'cython',
       'setuptools_scm',
    ],
    scripts=['scripts/utm-converter'],
    test_suite='test',
)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
