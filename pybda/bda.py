# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import
import time
import numpy
from oskar.bda import BDA


def run(config, vis, input_amp_name):
    num_times = config['sim']['observation']['num_times']
    time_inc = config['sim']['observation']['dump_time_s']
    freq_hz = config['sim']['observation']['freq_hz']
    max_fact = config['baseline_average']['max_fact']
    fov_radius_deg = config['baseline_average']['fov_radius_deg']
    max_average_time_s = config['baseline_average']['max_average_time_s']
    num_antennas = vis['num_antennas']
    num_baselines = vis['num_baselines']
    num_input_vis = len(vis['uu'])
    print('- Applying baseline-dependent time averaging...')
    print('  - Max factor              : %.4f' % max_fact)
    print('  - FoV radius              : %.2f deg' % fov_radius_deg)
    print('  - Max average time        : %.1f s' % max_average_time_s)
    print('  - No. input visibilities  : %i' % num_input_vis)

    # Create and set up the averager.
    t0 = time.time()
    h = BDA(num_antennas)
    duvw_max = h.set_compression(max_fact, fov_radius_deg, 299792458.0/freq_hz, 
        max_average_time_s)
    print('  - Delta UVW max           : %.4f m' % duvw_max)
    h.set_delta_t(time_inc)
    h.set_num_times(num_times)
    h.set_initial_coords(vis['uu'][0:num_baselines],
        vis['vv'][0:num_baselines], vis['ww'][0:num_baselines])
    for t in range(num_times):
        uu = None
        vv = None
        ww = None
        i0 = t * num_baselines
        i1 = (t + 1) * num_baselines
        i2 = (t + 2) * num_baselines
        amp = vis[input_amp_name][i0:i1]
        if t < num_times - 1:
            uu = vis['uu'][i1:i2]
            vv = vis['vv'][i1:i2]
            ww = vis['ww'][i1:i2]
        h.add_data(t, amp, uu, vv, ww)
    ave_data = h.finalise()
    num_output_vis = len(ave_data['uu'])
    compression_ratio = float(num_input_vis) / float(num_output_vis)
    print('  - No. output visibilities : %i' % num_output_vis)
    print('  - Compression ratio       : %.3f' % compression_ratio)
    print('  - Visibilities averaged in %.2f s' % (time.time() - t0))

    return ave_data, compression_ratio
