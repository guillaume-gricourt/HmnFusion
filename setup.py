#!/usr/bin/env python
"""Analysis fusion from DNA genomics"""
import os
import sys

from setuptools import setup, find_packages

setup_args = {}

install_requires=[
	'cnvkit >= 0.9.6',
	'django >= 1.11.25',
	'Jinja2 >= 2.10.3',
	'networkx >= 2.4.0',
	'seaborn >= 0.9.0',
	'pdfkit >= 0.6.1',
	'pysam >= 0.10.0',
	'xlsxwriter >= 1.1.6',
]

DIR = (os.path.dirname(__file__) or '.')
VERSION = ""
with open(os.path.join(DIR, "hmnfusion", '_version.py')) as fid:
	VERSION = fid.readline().split("=")[-1].strip().replace('"','')

setup_args.update(
	name='HmnFusion',
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
	packages=find_packages(),
	scripts = [os.path.join(DIR, "HmnFusion")],
	install_requires=install_requires,
	python_requires='>=3',
)
setup(**setup_args)

