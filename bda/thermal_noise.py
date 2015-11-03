# -*- coding: utf-8 -*-

import math


# https://www.skatelescope.org/wp-content/uploads/2012/07/SKA-TEL-SKO-DD-001-1_BaselineDesign1.pdf
# Table 7, p49


k0 = 1.38064852e-23
c0 = 299792458.0  # m/s
nu = 700.0e6  # Hz
wavelength = c0 / nu  # m
dish_diameter = 15.0  # m
dish_area = math.pi * (dish_diameter / 2)**2
eta = 0.65

A_eff = dish_area
T_sys = 60.0 * wavelength**2.55
T_sys = 28.0
S_sys = 2.0 * k0 * T_sys / (A_eff * eta)
S_sys /= 1.0e-26

print 'SEFD:', S_sys, 'Jy'

# bw = 700.0e6 / 2**16
bw = 700.0e6 / 1
dt = 0.1

sigma_pq = S_sys / (2.0 * bw * dt)**0.5
sigma_pq /= 2.0**0.5

print 'sigma_pq =', sigma_pq, 'Jy'

num_antennas = 133
num_baselines = (num_antennas * (num_antennas - 1)) / 2
num_times = 10
sigma_im = sigma_pq / (0.9*(num_baselines * num_times)**0.5)

print num_baselines
print 'sigma_im =', sigma_im, 'Jy/beam'
print 'sigma_im =', sigma_im * 1e6, 'uJy/beam'

# print S_sys
# print (674.0 / (1.0 * (2.0 * bw * num_baselines)**0.5)) * 1e3
