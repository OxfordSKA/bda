# -*- coding: utf-8 -*-
"""Prototype of python visibility simulator"""

from __future__ import print_function, absolute_import
from pyuvwsim import (load_station_coords, convert_enu_to_ecef,
                      evaluate_baseline_uvw)
from os.path import join
import math
import numpy
import time
import matplotlib.pyplot as pyplot


if __name__ == '__main__':
    freq_hz = 700.0e6
    lon = math.radians(21.442909)
    lat = math.radians(-30.739475)
    alt = 0.0
    source_ra = math.radians(-90.35458487600000)
    source_dec = math.radians(-7.67112399060000)
    ra = math.radians(-90.3545848760)
    dec = math.radians(-8.5711239906)
    mjd_mid = 57086.113194
    obs_length_s = 5.0
    dump_time_s = 0.1
    over_sample = 10
    dt_s = dump_time_s / over_sample
    num_times = int(obs_length_s / dt_s)
    print('- num times', num_times)

    # TODO convert source ra, dec to l, m, n
    l = math.sin(math.radians(1.0))
    m = 0.0
    n = (1 - l**2 - m**2)**0.5

    telescope_model = join('../models', 'ska1_meerkat_mid_combined_july_2015.tm')
    x, y, z = load_station_coords(join(telescope_model, 'layout.txt'))
    x, y, z = convert_enu_to_ecef(x, y, z, lon, lat, alt)

    num_antennas = x.shape[0]
    num_baselines = num_antennas * (num_antennas - 1) / 2
    num_coords = num_times * num_baselines
    print('- num antennas', num_antennas)
    print('- num baselines', num_baselines)
    print('- num coords', num_coords)
    print('- coord memory = %.2f MB'
           % ((num_coords * 3 * 8) / 1024.0**2))

    uu = numpy.zeros(num_coords, dtype='f8')
    vv = numpy.zeros(num_coords, dtype='f8')
    ww = numpy.zeros(num_coords, dtype='f8')
    t0 = time.time()
    for i in range(num_times):
        t = -obs_length_s / 2.0 + i * dt_s + dt_s / 2.0
        mjd = mjd_mid + (t / 86400.0)
        i0 = i * num_baselines
        i1 = i0 + num_baselines
        uu_, vv_, ww_ = evaluate_baseline_uvw(x, y, z, ra, dec, mjd)
        uu[i0:i1] = uu_
        vv[i0:i1] = vv_
        ww[i0:i1] = ww_
    print('- Time taken to generate coords = %.3fs' % (time.time() - t0))

    t0 = time.time()
    wavelength = 299792458.0 / freq_hz
    k = 2.0 * math.pi / wavelength
    amp = numpy.exp(1.0j * k * (uu * l + vv * m + ww * n))
    print('- Time taken to generate amplitudes = %.3fs' % (time.time() - t0))

    uv_dist = (uu**2 + vv**2)**0.5
    x = uv_dist[::100]
    y = amp[::100].real

    print('- plotting...')
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, '+')
    pyplot.show()




