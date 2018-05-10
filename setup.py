"""Poropyck setup file"""
from setuptools import setup, Extension
import numpy

EXT_MODULES = [Extension(
    name='poropyck.dtw',
    sources=['poropyck/src/dtw.pyx'],
    include_dirs=[numpy.get_include()],
    extra_compile_args=['-O3'],
    language='c')]

setup(
    name='poropyck',
    version='1.2.0',
    author='Evert Duran Quintero',
    author_email='edur409@aucklanduni.ac.nz',
    packages=['poropyck'],
    include_package_data=True,
    ext_modules=EXT_MODULES,
    python_requires='~=3.4',
    install_requires=[
        'numpy>=1.13.3',
        'matplotlib>=2.1.0',
        'scipy>=0.19.0',
        'mcerp3>=1.0.2'],
    entry_points={'console_scripts': [
        'install_poropyck_samples = poropyck.install_samples:install_samples',
        'simple_harmonics_1 = poropyck.simple_harmonics:simple_harmonics_1',
        'simple_harmonics_2 = poropyck.simple_harmonics:simple_harmonics_2',
        'pick_vp_dtw = poropyck.Pick_Vp_DTW:pick_vp_dtw',
        'pick_vs_dtw = poropyck.Pick_Vs_DTW:pick_vs_dtw',
        'cross_corr_p = poropyck.Cross_corr_P:cross_corr_p',
        'cross_corr_s = poropyck.Cross_corr_S:cross_corr_s',
        'vp_dtw_saturated_loop = poropyck.Vp_DTW_saturated_loop:vp_dtw_saturated_loop',
        'vs_dtw_saturated_loop = poropyck.Vs_DTW_saturated_loop:vs_dtw_saturated_loop',
        'densities_porosity = poropyck.Densities_Porosity:densities_porosity',
        'densities_porosity2 = poropyck.Densities_Porosity2:densities_porosity2',
        'elastic_constants = poropyck.Elastic_Constants:elastic_constants',
        'elastic_constants_pyc = poropyck.Elastic_Constants_Pyc:elastic_constants',
        'pycnometer_densities = poropyck.Pycnometer_densities:pycnometer_densities'],},
    )
