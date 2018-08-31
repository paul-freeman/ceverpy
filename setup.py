"""Poropyck setup file"""
from setuptools import setup

setup(
    name='poropyck',
    version='1.3.2',
    author='Evert Duran Quintero',
    author_email='edur409@aucklanduni.ac.nz',
    packages=['poropyck'],
    include_package_data=True,
    entry_points={'console_scripts': [
        'pick_dtw = poropyck.pick_dtw:pick_velocity'], },
)
