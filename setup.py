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
                  sources=["cy_src/sam_functions.pyx"])
    ]
    cmd_class.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("pyMiniSam",
                  sources=["cy_src/sam_functions.pyx"])
    ]

setup(name='pyMiniSam',
      setup_requires=['cython'],
      version=1.0,
      author="apfejes",
      long_description=open('README.txt').read(),
      cmdclass={'build_ext': build_ext},
      packages=find_packages(),
      zip_safe=False,
      include_package_data=True,
      ext_modules=cythonize(ext_modules, annotate=True))
