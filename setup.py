from distutils.core import setup
from Cython.Build import cythonize

setup(name='pyMiniSam',
      ext_modules=cythonize("sam_functions.pyx", annotate=True))
