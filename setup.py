#!/usr/bin/env python

from setuptools import setup, find_packages
from os import path

setup(
    name="spgr",
    version="0.1.0",
    description="Simple python space group library",

    author="Stef Smeets",
    author_email="stef.smeets@mat.ethz.ch",
    license="GPL",
    url="https://github.com/stefsmeets/spgr",

    classifiers=[
        'Programming Language :: Python :: 2.7',
    ],

    packages=["spgr"],

    install_requires=["numpy"],

    package_data={
        "": ["LICENCE", "readme.md"],
        "spgr": ["*.py"]
    },

    entry_points={
        'console_scripts': [
            'spgr = spgr.app:main',
        ]
    }

)
