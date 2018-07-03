from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

USE_CYTHON = True

cmd_class = {}
ext_modules = []

if USE_CYTHON:
    try:
        from Cython.Distutils import build_ext
    except ImportError:
        if USE_CYTHON == 'auto':
            USE_CYTHON = False
        else:
            raise

if USE_CYTHON:
    ext_modules += [
        Extension("pyMiniSam.sam_functions",
                  sources=["sam_functions.pyx"],
                  librarires=),
    ]
    cmd_class.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("pyMiniSam.sam_functions",
                  sources=["sam_functions.pyx"],
                  librarires=),
    ]

setup(name='pyMiniSam',
      version=1.0,
      author="apfejes",
      long_description=open('README.txt').read(),
      ext_modules=cythonize(ext_modules, annotate=True))
