from distutils.core import setup
from Cython.Build import cythonize

setup(name='pyMiniSam',
      version=1.0,
      author="apfejes",

      ext_modules=cythonize("sam_functions.pyx", annotate=True))
