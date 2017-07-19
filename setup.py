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
  name='GooseSolid',
  version=__version__,
  author='Tom de Geus',
  author_email='tom@geus.me',
  url='https://github.com/tdegeus/GooseSolid',
  description='Material models',
  long_description='',
  license='MIT',
  ext_modules=ext_modules,
  install_requires=['pybind11>=2.1.0','cppmat>=0.1.7','goosempl>=0.1.3'],
  cmdclass={'build_ext': cppmat.BuildExt},
  zip_safe=False,
)
