#!/usr/bin/env python
"""Analysis fusion from DNA genomics"""
import os
import sys

from setuptools import setup, find_packages

setup_args = {}

install_requires=[
    'beautifulsoup4>=4.9.0',
    'natsort>=7.0.0',
    'numpy>=1.18.3',
    'pandas>=1.0.3',
    'pysam>=0.15.4',
    'et-xmlfile',
    'openpyxl'
]

# Cmd in meta.yml
#script: "{{ PYTHON }} -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv"
    
DIR = (os.path.dirname(__file__) or '.')
NAME, VERSION = '', ''
with open(os.path.join(DIR, 'hmnfusion', '_version.py')) as fid:
    lines = fid.read().splitlines()
    NAME = lines[0].split("=")[-1].strip().replace('"','')
    VERSION = lines[1].split("=")[-1].strip().replace('"','')

setup_args.update(
    name=NAME,
    version=VERSION,
    description=__doc__,
    author='Guillaume Gricourt',
    author_email='guillaume.gricourt@aphp.fr',
    #url='',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    packages=['hmnfusion'],
    entry_points = {
        'console_scripts' : [ 'hmnfusion=hmnfusion.__main__:main']
    },
    install_requires=install_requires,
    python_requires='>=3.5,<=3.6.13',
)
setup(**setup_args)

