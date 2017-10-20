desc = '''
GooseMaterial is a C++ module, wrapped in Python, that provides several constitutive models for
solids. It can be used for example in finite element programs (at the integration point level).
'''

from setuptools import setup, Extension

import sys,re
import setuptools
import pybind11
import cppmat

header = open('../src/GooseMaterial/Macros.h','r').read()
world  = re.split('(.*)(\#define GOOSEMATERIAL_WORLD_VERSION\ )([0-9]+)(.*)',header)[3]
major  = re.split('(.*)(\#define GOOSEMATERIAL_MAJOR_VERSION\ )([0-9]+)(.*)',header)[3]
minor  = re.split('(.*)(\#define GOOSEMATERIAL_MINOR_VERSION\ )([0-9]+)(.*)',header)[3]

__version__ = '.'.join([world,major,minor])

ext_modules = [
  Extension(
    'GooseMaterial',
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
  name               = 'GooseMaterial',
  description        = 'Material models in C++ (and Python)',
  long_description   = desc,
  keywords           = 'FEM, Mechanics, C++, C++11, Python bindings, pybind11',
  version            = __version__,
  license            = 'GPLv3',
  author             = 'Tom de Geus',
  author_email       = 'tom@geus.me',
  url                = 'https://github.com/tdegeus/GooseMaterial',
  ext_modules        = ext_modules,
  extra_compile_args = ["-DNDEBUG"], # switch off assertions
  install_requires   = ['pybind11>=2.2.0','cppmat>=0.2.14','goosempl>=0.1.3'],
  cmdclass           = {'build_ext': cppmat.BuildExt},
  zip_safe           = False,
)
