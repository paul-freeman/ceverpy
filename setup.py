"""ceverpy setup file"""
from setuptools import setup

setup(
    name='ceverpy',
    version='1.0.0',
    author='Evert Duran Quintero',
    author_email='edur409@aucklanduni.ac.nz',
    packages=['ceverpy'],
    install_requires=[
        'numpy>=1.13.3',
        'matplotlib>=2.1.0',
        'scipy>=0.19.0',
        'mcerp3>=1.0.0',
        'rpy2>=2.9.0'],
    entry_points={'console_scripts': [
        'simple_harmonics = ceverpy.Simple_harmonics:simple_harmonics',
        'pick_vp_dtw = ceverpy.Pick_Vp_DTW:pick_vp_dtw',
        'pick_vs_dtw = ceverpy.Pick_Vs_DTW:pick_vs_dtw',
        'cross_corr_p = ceverpy.Cross_corr_P:cross_corr_p',
        'cross_corr_s = ceverpy.Cross_corr_S:cross_corr_s',
        'vp_dtw_saturated_loop = ceverpy.Vp_DTW_saturated_loop:vp_dtw_saturated_loop',
        'vs_dtw_saturated_loop = ceverpy.Vs_DTW_saturated_loop:vs_dtw_saturated_loop',
        'densities_porosity = ceverpy.Densities_Porosity:densities_porosity',
        'densities_porosity2 = ceverpy.Densities_Porosity2:densities_porosity2',
        'elastic_constants = ceverpy.Elastic_Constants:elastic_constants',
        'elastic_constants_pyc = ceverpy.Elastic_Constants_Pyc:elastic_constants_pyc'],},
    )
