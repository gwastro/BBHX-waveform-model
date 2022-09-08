#!/usr/bin/env python
"""
setup.py file for BBHX waveform into pycbc waveform plugin package
"""

from setuptools import Extension, setup, Command
from setuptools import find_packages

VERSION = '0.0.dev0'

setup (
    name = 'pycbc-bbhx',
    version = VERSION,
    description = 'Plugin of BBHX fast TDI waveform generator into PyCBC',
    long_description = open('descr.rst').read(),
    author = 'Connor Weaving',
    author_email = 'connor.weaving@port.ac.uk',
    url = 'http://www.pycbc.org/',
    download_url = 'https://github.com/gwastro/revchirp/tarball/v%s' % VERSION,
    keywords = ['pycbc', 'signal processing', 'gravitational waves'],
    install_requires = ['pycbc'],
    py_modules = ['bbhx'],
    entry_points = {"pycbc.waveform.fd":"bbhx = bbhx:bbhx_fd"},
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
