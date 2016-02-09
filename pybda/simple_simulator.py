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
import oskar


def convert_ra_dec_to_relative_directions(ra, dec, ra0, dec0):
    d_ra = ra - ra0
    cos_d_ra = numpy.cos(d_ra)
    sin_d_ra = numpy.sin(d_ra)
    cos_dec  = numpy.cos(dec)
    sin_dec  = numpy.sin(dec)
    cos_dec0 = math.cos(dec0)
    sin_dec0 = math.sin(dec0)
    l = cos_dec * sin_d_ra
    m = cos_dec0 * sin_dec - sin_dec0 * cos_dec * cos_d_ra
    n = sin_dec0 * sin_dec + cos_dec0 * cos_dec * cos_d_ra
    return l, m, n


def generate_baseline_uvw(x, y, z, ra_rad, dec_rad, \
            freq_hz, mjd_mid, length_s, dt_s):
    uu = numpy.zeros(num_coords, dtype='f8')
    vv = numpy.zeros(num_coords, dtype='f8')
    ww = numpy.zeros(num_coords, dtype='f8')
    for i in range(num_times):
        t = -length_s / 2.0 + i * dt_s + dt_s / 2.0
        mjd = mjd_mid + (t / 86400.0)
        i0 = i * num_baselines
        i1 = i0 + num_baselines
        uu_, vv_, ww_ = evaluate_baseline_uvw(x, y, z, ra_rad, dec_rad, mjd)
        uu[i0:i1] = uu_
        vv[i0:i1] = vv_
        ww[i0:i1] = ww_
    wavelength = 299792458.0 / freq_hz
    uu /= wavelength
    vv /= wavelength
    ww /= wavelength
    return uu, vv, ww


if __name__ == '__main__':
    freq_hz = 700.0e6
    lon = math.radians(21.442909)
    lat = math.radians(-30.739475)
    alt = 0.0
    source_ra = math.radians(-90.35458487600000)
    source_dec = math.radians(-7.67112399060000)
    vis_ra = math.radians(-90.3545848760)
    vis_dec = math.radians(-8.5711239906)
    img_ra = math.radians(-90.35458487600000)
    img_dec = math.radians(-7.67112399060000)
    mjd_mid = 57086.113194
    obs_length_s = 5.0
    dump_time_s = 0.1
    over_sample = 10
    dt_s = dump_time_s / over_sample
    num_times = int(obs_length_s / dt_s)
    print('- num times', num_times)

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

    # Generate actual UVW coordinates in wavelengths.
    t0 = time.time()
    uu, vv, ww = generate_baseline_uvw(x, y, z, vis_ra, vis_dec, \
            freq_hz, mjd_mid, obs_length_s, dt_s)
    print('- Time taken to generate coords = %.3fs' % (time.time() - t0))

    # Generate amplitudes.
    t0 = time.time()
    l, m, n = convert_ra_dec_to_relative_directions(source_ra, source_dec, \
        vis_ra, vis_dec)
    amp = numpy.exp(1.0j * 2.0 * math.pi * (uu * l + vv * m + ww * (n - 1.0)))
    print('- Time taken to generate amplitudes = %.3fs' % (time.time() - t0))

    # Phase rotate amplitudes.
    t0 = time.time()
    l, m, n = convert_ra_dec_to_relative_directions(img_ra, img_dec, 
        vis_ra, vis_dec)
    delta_l = 0.0 - l
    delta_m = 0.0 - m
    delta_n = 1.0 - n
    amp *= numpy.exp(1.0j * 2.0 * math.pi * \
        (uu * delta_l + vv * delta_m + ww * delta_n))
    print('- Time taken to phase rotate amplitudes = %.3fs' % (time.time() - t0))

    # Generate rotated UVW coordinates in wavelengths.
    t0 = time.time()
    uu, vv, ww = generate_baseline_uvw(x, y, z, img_ra, img_dec, \
            freq_hz, mjd_mid, obs_length_s, dt_s)
    print('- Time taken to generate rotated coords = %.3fs' % (time.time() - t0))

    # Make image.
    t0 = time.time()
    im = oskar.image.make(uu, vv, ww, amp, 0.1, 2048)
    print('- Time taken to generate image = %.3fs' % (time.time() - t0))
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.imshow(im, interpolation='nearest')
    pyplot.show()

    uv_dist = (uu**2 + vv**2)**0.5
    x = uv_dist[::100]
    y = amp[::100].real

    print('- plotting...')
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, '+')
    pyplot.show()




