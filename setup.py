#!/usr/bin/env python

from setuptools import setup, find_packages
from os import path

setup(
    name="xcore",
    version="0.4.0",
    description="Simple python crystallography library",

    author="Stef Smeets",
    author_email="stef.smeets@mmk.su.se",
    license="GPL",
    url="https://github.com/stefsmeets/xcore",

    classifiers=[
        'Programming Language :: Python :: 2.7',
    ],

    packages=["xcore", "xcore.scattering", "xcore.spacegroup_tables"],

    install_requires=["numpy", "pandas"],

    package_data={
        "": ["LICENCE", "readme.md"],
        "xcore": ["*.py", "spacegroups.txt"]
    },

    entry_points={
        'console_scripts': [
            'spgr    = xcore.app:main',
            'cif2hkl = xcore.formats:cif2hkl_entry',
            'make_focus = xcore.formats:make_focus_entry'
        ]
    }

)
