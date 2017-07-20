#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
import os, sys

srclist = ['sgglobal.c','sgcb.c','sgcharmx.c','sgfile.c',
           'sggen.c','sghall.c','sghkl.c','sgltr.c','sgmath.c','sgmetric.c',
           'sgnorm.c','sgprop.c','sgss.c','sgstr.c','sgsymbols.c',
           'sgtidy.c','sgtype.c','sgutil.c','runtests.c','sglitemodule.c']

srclist = [os.path.join('src/sglite', f) for f in srclist]
sglite_ext = Extension(
    'xcore.sglite',
    sources=srclist,
    define_macros=[('PythonTypes', 1)]
    )

setup(
    name="xcore",
    version="0.6.0",
    description="Simple python crystallography library",

    author="Stef Smeets",
    author_email="stef.smeets@mmk.su.se",
    license="GPL",
    url="https://github.com/stefsmeets/xcore",

    classifiers=[
        'Programming Language :: Python :: 2.7',
    ],

    packages=["xcore", "xcore.scattering"],

    ext_modules = [sglite_ext],

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
            'make_superflip = xcore.formats:make_superflip_entry'
        ]
    }

)
