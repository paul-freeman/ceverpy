"""Created on Tue Feb 6 15:44:40 2018

@author: leo
"""
from .poropyck_common import simple_harmonics


def simple_harmonics_1():
    """simple harmonics - version 1"""
    simple_harmonics(a2_value=[1.0, 1.5, 2.0, 3.0, 1.0, 0.0, 0.0],
                     cross_correlation=True,
                     use_manager=False,
                     print_time_lag=True,
                     labels=True)


def simple_harmonics_2():
    """simple harmonics - version 2"""
    simple_harmonics(a2_value=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     cross_correlation=False,
                     use_manager=True,
                     print_time_lag=False,
                     labels=False)
