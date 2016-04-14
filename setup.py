#!/usr/bin/env python

from setuptools import setup, find_packages
from os import path

setup(
    name="spgr",
    version="0.2.0",
    description="Simple python crystallography library",

    author="Stef Smeets",
    author_email="stef.smeets@mmk.su.se",
    license="GPL",
    url="https://github.com/stefsmeets/spgr",

    classifiers=[
        'Programming Language :: Python :: 2.7',
    ],

    packages=["spgr", "spgr.scattering", "spgr.spacegroup"],

    install_requires=["numpy", "pandas"],

    package_data={
        "": ["LICENCE", "readme.md"],
        "spgr": ["*.py", "spacegroups.txt"]
    },

    entry_points={
        'console_scripts': [
            'spgr    = spgr.app:main',
            'cif2hkl = spgr.diffraction.cif2hkl_entry'
        ]
    }

)
