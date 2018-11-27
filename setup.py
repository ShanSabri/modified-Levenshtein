from distutils.core import setup
from Cython.Build import cythonize
# import numpy

setup(
        ext_modules=cythonize("cLev.pyx")
        # include_dirs=[numpy.get_include()] # needed for np support
)
