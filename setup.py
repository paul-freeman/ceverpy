"""Poropyck setup file"""
from setuptools import setup

setup(
    name='poropyck',
    version='1.3.0',
    author='Evert Duran Quintero',
    author_email='edur409@aucklanduni.ac.nz',
    packages=['poropyck'],
    include_package_data=True,
    python_requires='~=3.4',
    install_requires=[
        'numpy>=1.13.3,<1.15',
        'matplotlib>=2.1.0',
        'scipy>=0.19.0',
        'dtw>=1.0.0',
        'mcerp3>=1.0.2'],
    entry_points={'console_scripts': [
        'poropyck = poropyck.pick_dtw:main'], },
)
