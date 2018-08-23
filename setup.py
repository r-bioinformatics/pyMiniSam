from setuptools import setup, find_packages
from Cython.Build import cythonize
from setuptools.extension import Extension
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
                  sources=["pyMiniSam/sam_functions.pyx"])
    ]
    cmd_class.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("pyMiniSam",
                  sources=["pyMiniSam/sam_functions.pyx"])
    ]

setup(name='pyMiniSam',
      setup_requires=['cython'],
      version=1.0,
      author="apfejes",
      long_description=open('README.txt').read(),
      cmdclass={'build_ext': build_ext},
      packages=['pyMiniSam'],
      zip_safe=False,
      ext_modules=cythonize(ext_modules, annotate=True),
      setup_requires = ['Cython'],
      test_deps = ['Cython'],
)
