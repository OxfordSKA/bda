# -*- coding: utf-8 -*-
"""Simple visibility simulator in python."""

from __future__ import print_function, absolute_import
from pyuvwsim import (load_station_coords, convert_enu_to_ecef,
                      evaluate_baseline_uvw)
from os.path import join
import math
import numpy
from psutil import virtual_memory
from pybda.corrupting_gains import eval_complex_gains
import time
from oskar._bda_utils import apply_gains, check_ref_count


def convert_ra_dec_to_relative_directions(ra, dec, ra0, dec0):
    d_ra = ra - ra0
    cos_d_ra = numpy.cos(d_ra)
    sin_d_ra = numpy.sin(d_ra)
    cos_dec = numpy.cos(dec)
    sin_dec = numpy.sin(dec)
    cos_dec0 = math.cos(dec0)
    sin_dec0 = math.sin(dec0)
    l = cos_dec * sin_d_ra
    m = cos_dec0 * sin_dec - sin_dec0 * cos_dec * cos_d_ra
    n = sin_dec0 * sin_dec + cos_dec0 * cos_dec * cos_d_ra
    return l, m, n


def generate_baseline_uvw(x, y, z, ra_rad, dec_rad, num_times, num_baselines,
                          mjd_start, dt_s):
    num_coords = num_times * num_baselines
    uu = numpy.zeros(num_coords, dtype='f8')
    vv = numpy.zeros(num_coords, dtype='f8')
    ww = numpy.zeros(num_coords, dtype='f8')
    for i in range(num_times):
        t = i * dt_s + dt_s / 2.0
        mjd = mjd_start + (t / 86400.0)
        i0 = i * num_baselines
        i1 = i0 + num_baselines
        uu_, vv_, ww_ = evaluate_baseline_uvw(x, y, z, ra_rad, dec_rad, mjd)
        uu[i0:i1] = uu_
        vv[i0:i1] = vv_
        ww[i0:i1] = ww_
    return uu, vv, ww


def simulate(config):
    telescope = config['sim']['telescope']
    obs = config['sim']['observation']
    freq_hz = obs['freq_hz']
    lon = math.radians(telescope['lon_deg'])
    lat = math.radians(telescope['lat_deg'])
    alt = telescope['alt_m']
    vis_ra = math.radians(obs['ra_deg'])
    vis_dec = math.radians(obs['dec_deg'])
    mjd_start = obs['start_time_mjd']
    num_dumps = obs['num_times']
    dump_time_s = obs['dump_time_s']
    over_sample = obs['over_sample']
    dt_s = dump_time_s / over_sample
    num_times = num_dumps * over_sample
    telescope_model = telescope['path']
    x, y, z = load_station_coords(join(telescope_model, 'layout.txt'))
    x, y, z = convert_enu_to_ecef(x, y, z, lon, lat, alt)
    num_antennas = x.shape[0]
    num_baselines = num_antennas * (num_antennas - 1) / 2
    num_vis = num_times * num_baselines
    sky_model = config['sim']['sky_file']
    sky = numpy.loadtxt(sky_model, delimiter=',')
    sky = sky.reshape((-1, 3))
    l, m, n = convert_ra_dec_to_relative_directions(sky[:, 0], sky[:, 1],
                                                    vis_ra, vis_dec)
    amp = numpy.zeros(num_vis, dtype='c16')

    print('- Simple simulator.')
    print('  - No. antennas  =', num_antennas)
    print('  - No. baselines =', num_baselines)
    print('  - Obs. length   = %.1f s' % (num_dumps * dump_time_s))
    print('  - No. times     = %i (no. dumps: %i, over-sample: %i)' %
          (num_times, num_dumps, over_sample))
    print('  - No. vis       =', num_vis)
    print('  - No. sources   =', sky.shape[0])

    num_bytes = num_vis * 7 * 8
    mem = virtual_memory()
    print('  - Memory required = %.1f / %.1f MB' %
          (num_bytes / 1024.0**2, mem.total / 1024.0**2))
    if num_bytes >= mem.total:
        raise RuntimeError('Not enough system memory for requested '
                           'simulation.')

    # Generate UVW coordinates.
    t0 = time.time()
    uu, vv, ww = generate_baseline_uvw(x, y, z, vis_ra, vis_dec, num_times,
                                       num_baselines, mjd_start, dt_s)
    t_coords = time.time() - t0

    # Generate amplitudes.
    t1 = time.time()
    wavelength = 299792458.0 / freq_hz
    for i in range(len(l)):
        phase = (2.0 * math.pi / wavelength) * \
                (uu * l[i] + vv * m[i] + ww * (n[i] - 1.0))
        amp += numpy.exp(1.0j * phase)
    t_amp = time.time() - t1
    print('  - Total simulation time = %.2f s (coords: %.2f s, amp: %.2f s)'
          % (time.time() - t0, t_coords, t_amp))

    return {'model': amp, 'uu': uu, 'vv': vv, 'ww': ww}


def simulate_2(config):
    """Simulate with and without corruptions followed by averaging.
    Simulation to be performed in blocks due to memory constraints.
    """
    telescope = config['sim']['telescope']
    obs = config['sim']['observation']
    corrupt = config['corrupt']
    freq_hz = obs['freq_hz']
    wavelength = 299792458.0 / freq_hz
    lon = math.radians(telescope['lon_deg'])
    lat = math.radians(telescope['lat_deg'])
    alt = telescope['alt_m']
    ra = math.radians(obs['ra_deg'])
    dec = math.radians(obs['dec_deg'])
    mjd_start = obs['start_time_mjd']
    num_dumps = obs['num_times']
    dump_time_s = obs['dump_time_s']
    over_sample = obs['over_sample']
    dt_s = dump_time_s / over_sample
    num_times = num_dumps * over_sample
    telescope_model = telescope['path']
    x, y, z = load_station_coords(join(telescope_model, 'layout.txt'))
    x, y, z = convert_enu_to_ecef(x, y, z, lon, lat, alt)
    num_antennas = x.shape[0]
    num_baselines = num_antennas * (num_antennas - 1) / 2
    num_vis = num_dumps * num_baselines
    sky_model = config['sim']['sky_file']
    sky = numpy.loadtxt(sky_model, delimiter=',')
    sky = sky.reshape((-1, 3))
    num_sources = sky.shape[0]
    source_ra = numpy.radians(sky[:, 0])
    source_dec = numpy.radians(sky[:, 1])

    l, m, n = convert_ra_dec_to_relative_directions(
        source_ra, source_dec, ra, dec)
    tau = corrupt['tau_s']
    hurst_amp = corrupt['amplitude']['hurst']
    adev_amp = corrupt['amplitude']['allan_dev']
    std_t_mid_amp = corrupt['amplitude']['std_t_mid']
    hurst_phase = corrupt['phase']['hurst']
    adev_phase = corrupt['phase']['allan_dev']
    std_t_mid_phase = corrupt['phase']['std_t_mid']
    smoothing_length = corrupt['smoothing_length']

    print('- Simple simulator.')
    print('  - No. antennas  =', num_antennas)
    print('  - No. baselines =', num_baselines)
    print('  - Obs. length   = %.1f s' % (num_dumps * dump_time_s))
    print('  - No. times     = %i (no. dumps: %i, over-sample: %i)' %
          (num_times, num_dumps, over_sample))
    print('  - No. vis       =', num_vis)
    print('  - No. sources   =', num_sources)
    print('  - Corruptions:')
    print('    -  Hurst amp %.1f, phase %.1f' % (hurst_amp, hurst_phase))
    print('    -  A. dev amp %.1e, phase %.1e' % (adev_amp, adev_phase))

    num_bytes = num_vis * 8 * 7 + (num_antennas * num_times) * 16
    mem = virtual_memory()
    print('  - Mem. required = %.1f / %.1f MB' %
          (num_bytes / 1024.0**2, mem.total / 1024.0**2))
    if num_bytes >= mem.total:
        raise RuntimeError('Not enough system memory for requested '
                           'simulation.')

    # Generate corruptions
    gains = numpy.empty((num_antennas, num_times), dtype='c16')
    t0 = time.time()
    for i in range(num_antennas):
        gains[i, :] = eval_complex_gains(num_times, dt_s, hurst_amp, adev_amp,
                                         std_t_mid_amp, hurst_phase, adev_phase,
                                         std_t_mid_phase, smoothing_length,
                                         tau)
    conj_gains = numpy.conj(gains)
    print('  - Gains generated in %.1f s' % (time.time() - t0))

    # TODO-BM plot the gains ...

    # Simulation
    phase0 = 2.0 * math.pi / wavelength
    model = numpy.empty(num_vis, dtype='c16')
    data = numpy.empty(num_vis, dtype='c16')
    block_model = numpy.empty((over_sample, num_baselines), dtype='c16')
    block_data = numpy.empty((over_sample, num_baselines), dtype='c16')
    t1 = time.time()
    # Loop over correlator dumps
    for i in range(num_dumps):
        block_model.fill(0.0 + 0.0j)
        block_data.fill(0.0 + 0.0j)
        # Loop over times in a dump.
        for t in range(over_sample):
            delta_t = (i * dump_time_s) + (t * dt_s) + (dt_s / 2.0)
            mjd = mjd_start + (delta_t / 86400.0)
            uu_, vv_, ww_ = evaluate_baseline_uvw(x, y, z, ra, dec, mjd)
            for s in range(num_sources):
                phase = phase0 * (uu_ * l[s] + vv_ * m[s] + ww_ * (n[s] - 1.0))
                block_model[t, :] += numpy.exp(1.0j * phase)
            ig = i * over_sample + t
            block_data[t, :] = apply_gains(block_model[t, :], gains[:, ig])
            # data_ = apply_gains(block_model[t, :], gains[:, ig])
            # idx = 0
            # for p in range(num_antennas):
            #     for q in range(p + 1, num_antennas):
            #         gp = gains[p, ig]
            #         gq = conj_gains[q, ig]
            #         block_data[t, idx] = gp * block_model[t, idx] * gq
            #         idx += 1
            # if i < 3:
            #     print('t[%02i]' % t)
            #     print(' C:', data_[0])
            #     print(' P:', block_data[t, 0])
            #     print(' D:', numpy.max(numpy.abs(data_[:] - block_data[t, :])))
            # assert numpy.max(numpy.abs(data_[:] - block_data[t, :])) == 0.0

        # Average the block to get visibility data amplitudes for the dump.
        b0 = i * num_baselines
        b1 = b0 + num_baselines
        model[b0:b1] = numpy.mean(block_model, axis=0)
        data[b0:b1] = numpy.mean(block_data, axis=0)

    uu, vv, ww = generate_baseline_uvw(x, y, z, ra, dec, num_dumps,
                                       num_baselines, mjd_start, dump_time_s)

    print('  - Visibilities simulated in %.2f s' % (time.time() - t1))

    return {'model': model, 'data': data, 'uu': uu, 'vv': vv, 'ww': ww,
            'num_baselines': num_baselines, 'num_times': num_dumps,
            'num_antennas': num_antennas}


