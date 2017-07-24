#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
import os, sys
import glob

src_drc = "src\sglite"
headers = glob.glob(os.path.join(src_drc, "*.h"))

src_files = ['sgglobal.c','sgcb.c','sgcharmx.c','sgfile.c',
           'sggen.c','sghall.c','sghkl.c','sgltr.c','sgmath.c','sgmetric.c',
           'sgnorm.c','sgprop.c','sgss.c','sgstr.c','sgsymbols.c',
           'sgtidy.c','sgtype.c','sgutil.c','runtests.c','sglitemodule.c']

src_files = [os.path.join(src_drc, f) for f in src_files]
sglite_ext = Extension(
    'xcore.sglite',
    sources=src_files,
    define_macros=[('PythonTypes', 1)]
    )

setup(
    name="xcore",
    version="0.6.0",
    description="Crystallographic space group library in Python",

    author="Stef Smeets",
    author_email="stef.smeets@mmk.su.se",
    license="GPL",
    url="https://github.com/stefsmeets/xcore",

    classifiers=[
        'Programming Language :: Python :: 2.7',
    ],

    packages=["xcore", "xcore.scattering"],

    ext_modules = [sglite_ext],
    headers = headers,

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
