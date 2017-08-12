desc = '''
GooseSolid is a C++ module, wrapped in Python, that provides several constitutive models for solids.
It can be used for example in finite element programs (at the integration point level).
'''

from setuptools import setup, Extension

import sys
import setuptools
import pybind11
import cppmat

__version__ = '0.0.1'

ext_modules = [
  Extension(
    'GooseSolid',
    ['python_interface.cpp'],
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
      cppmat  .get_include(False),
      cppmat  .get_include(True )
    ],
    language='c++'
  ),
]

setup(
  name               = 'GooseSolid',
  description        = 'Material models in C++ (and Python)',
  long_description   = desc,
  keywords           = 'FEM, Mechanics, C++, C++11, Python bindings, pybind11',
  version            = __version__,
  license            = 'GPLv3',
  author             = 'Tom de Geus',
  author_email       = 'tom@geus.me',
  url                = 'https://github.com/tdegeus/GooseSolid',
  ext_modules        = ext_modules,
  extra_compile_args = ["-DNDEBUG"], # switch off assertions
  install_requires   = ['pybind11>=2.1.0','cppmat>=0.2.1','goosempl>=0.1.3'],
  cmdclass           = {'build_ext': cppmat.BuildExt},
  zip_safe           = False,
)
