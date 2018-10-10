from distutils.core import setup
from Cython.Build import cythonize
 
setup(
    ext_modules = cythonize("calibrate_model_cluster_per_set.pyx")
)