#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
doclink = """
Documentation
-------------

The full documentation is at http://rnamake_ens_gen.rtfd.org."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='rnamake_ens_gen',
    version='0.1.0',
    description='a series of functions to optimize a motif ensemble using rnamake and chemical mapping data',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    url='https://github.com/jyesselm/rnamake_ens_gen',
    packages=[
        'rnamake_ens_gen',
    ],
    package_dir={'rnamake_ens_gen': 'rnamake_ens_gen'},
    py_modules=[
        'rnamake_ens_gen / rnamake_ens_gen'
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='rnamake_ens_gen',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points = {
        'console_scripts' : [
        ]
    }
)
