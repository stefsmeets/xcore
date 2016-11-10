# from distutils.core import setup, Extension
from setuptools import Command, Extension, find_packages, setup

srclist = ['sgglobal.c','sgcb.c','sgcharmx.c','sgfile.c',
           'sggen.c','sghall.c','sghkl.c','sgltr.c','sgmath.c','sgmetric.c',
           'sgnorm.c','sgprop.c','sgss.c','sgstr.c','sgsymbols.c',
           'sgtidy.c','sgtype.c','sgutil.c','runtests.c','sglitemodule.c']

module = Extension('sglite', sources=srclist,
                   define_macros = [('PythonTypes', 1)])

setup (name='sglite',
       description = 'space group info',
       ext_modules = [module]
       )

