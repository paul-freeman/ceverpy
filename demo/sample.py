"""Example 1: basic poropyck script"""
import poropyck

lengths = [5.256, 5.25, 5.254, 5.254, 5.252, 5.252, 5.258, 5.265, 5.255, 5.252]
dry_data = 'NM11_2087_4A_dry.csv'
sat_data = 'NM11_2087_4A_sat.csv'

dtw = poropyck.DTW(dry_data, sat_data, lengths)
dry, sat = dtw.pick()

print('Dry:\n', dry)
print('Sat:\n', sat)
