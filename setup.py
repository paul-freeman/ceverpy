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
    version='1.2.1',
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
        'pick_vp_dtw = poropyck.pick_vp_dtw:pick_vp_dtw',
        'pick_vs_dtw = poropyck.pick_vs_dtw:pick_vs_dtw',
        'cross_corr_p = poropyck.cross_corr_p:cross_corr_p',
        'cross_corr_s = poropyck.cross_corr_s:cross_corr_s',
        'vp_dtw_saturated_loop = poropyck.vp_dtw_saturated_loop:vp_dtw_saturated_loop',
        'vs_dtw_saturated_loop = poropyck.vs_dtw_saturated_loop:vs_dtw_saturated_loop',
        'densities_porosity = poropyck.densities_porosity:densities_porosity',
        'densities_porosity2 = poropyck.densities_porosity2:densities_porosity2',
        'elastic_constants = poropyck.elastic_constants:elastic_constants',
        'elastic_constants_pyc = poropyck.elastic_constants_pyc:elastic_constants',
        'pycnometer_densities = poropyck.pycnometer_densities:pycnometer_densities'],},
    )
