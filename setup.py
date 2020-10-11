#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from distutils.core import setup
from setuptools import setup, find_packages

## conda skeleton can't find readme
import os
if os.path.isfile("README.MD"):
	with open("README.MD", "r") as fh:
		long_description = fh.read()
else:
	long_description="change-seq"

setup(
	name='changeseq',
	version='1.2.8', # fix ECDF importing error 
	description="Bioinformatic pipeline for the CHANGE-seq assay.",
	author="Shengdar Q Tsai, Martin Aryee, Ved V Topkar, Jose Malagon-Lopez",
	author_email='STSAI4@mgh.harvard.edu, Aryee.Martin@mgh.harvard.edu, vedtopkar@gmail.com, jose.lopez@mail.harvard.edu',
	url='https://github.com/tsailabSJ/changeseq',
	packages=['changeseq'],
	# package_dir={'changeseq':'changeseq'},
	license='LICENSE',
	scripts=['changeseq/changeseq.py','changeseq/alignReads.py','changeseq/visualization.py',
		'changeseq/callVariants.py','changeseq/findCleavageSites.py','changeseq/log.py',
		'changeseq/mergeReads.py','changeseq/referenceFree.py','changeseq/utility.py',
		'changeseq/validation.py'],
	package_data={'test': ["test/*"]},
	include_package_data=True,
	long_description=long_description,
	long_description_content_type='text/markdown'	,
	keywords='changeseq',
	classifiers=[
		'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Visualization',
		'Topic :: Scientific/Engineering :: Information Analysis',
		'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
		'Operating System :: Unix',
		'Natural Language :: English',
		"Programming Language :: Python :: 2",
		'Programming Language :: Python :: 2.6',
		'Programming Language :: Python :: 2.7'
	]
)
