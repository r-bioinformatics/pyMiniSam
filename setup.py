from setuptools import setup, find_packages
from Cython.Build import cythonize
from setuptools.extension import Extension
USE_CYTHON = True

def install_and_import(package):
    import importlib
    try:
        importlib.import_module(package)
    except ImportError:
        import pip
        pip.main(['install', package])
    finally:
        globals()[package] = importlib.import_module(package)


cmd_class = {}
ext_modules = []

if USE_CYTHON:
    try:
        install_and_import("cython")
        from Cython.Distutils import build_ext
    except ImportError or ModuleNotFoundError:
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
      install_requires=['Cython'],
)
