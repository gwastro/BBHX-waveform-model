#!/usr/bin/env python
"""
setup.py file for BBHX waveform into pycbc waveform plugin package
"""

from setuptools import Extension, setup, Command
from setuptools import find_packages

VERSION = '0.2'

setup (
    name = 'pycbc-bbhx-plugin',
    version = VERSION,
    description = 'Plugin of BBHX fast TDI waveform generator into PyCBC',
    author = ['Connor Weaving', 'Shichao Wu'],
    author_email = ['connor.weaving@port.ac.uk', 'shichao.wu@aei.mpg.de'],
    url = 'http://www.pycbc.org/',
    download_url = 'https://github.com/gwastro/BBHX-waveform-model/v%s' % VERSION,
    keywords = ['pycbc', 'signal processing', 'gravitational waves', 'lisa'],
    install_requires = ['pycbc'],
    py_modules = ['BBHX_Phenom'],
    entry_points = {"pycbc.waveform.fd_det":["BBHX_PhenomD=BBHX_Phenom:waveform_setup", "BBHX_PhenomHM=BBHX_Phenom:waveform_setup"],
                    "pycbc.waveform.fd_det_sequence":["BBHX_PhenomD=BBHX_Phenom:waveform_setup", "BBHX_PhenomHM=BBHX_Phenom:waveform_setup"],
                    "pycbc.waveform.length":"BBHX_PhenomD=BBHX_Phenom:imr_duration"},

    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
