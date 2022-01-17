#!/usr/bin/env python
"""Analysis fusion from DNA genomics"""
import glob
import yaml

from setuptools import setup


setup_args = {}

# Version
version = ''
fversion = glob.glob('**/_version.py')[0]
with open(fversion) as fid:
    lines = fid.read().splitlines()
    version = lines[0].split("=")[-1].strip().replace('"', '')

# App name - dependencies
env = {}
with open('environment.yml') as fid:
    env = yaml.safe_load(fid)
name = env['name']
install_requires = env['dependencies']

setup_args.update(
    name=name,
    version=version,
    author='Guillaume Gricourt',
    author_email='guipagui@gmail.com',
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
    packages=['hmnfusion'],
    entry_points={
        'console_scripts': ['hmnfusion=hmnfusion.__main__:main']
    },
    install_requires=install_requires,
)
setup(**setup_args)
