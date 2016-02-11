# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import
from pybda.simple_simulator import (convert_ra_dec_to_relative_directions,
                                    generate_baseline_uvw)
from pyuvwsim import (load_station_coords, convert_enu_to_ecef,
                      evaluate_baseline_uvw)
import time
import math
import numpy
from os.path import join
from oskar.image import make


def make_image(config, vis, image_id=0):
    telescope = config['sim']['telescope']
    image_config = config['imaging']['images'][image_id]
    obs_config = config['sim']['observation']
    mjd_start = obs_config['start_time_mjd']
    num_dumps = obs_config['num_times']
    dump_time_s = obs_config['dump_time_s']
    freq_hz = obs_config['freq_hz']
    vis_ra = math.radians(obs_config['ra_deg'])
    vis_dec = math.radians(obs_config['dec_deg'])
    image_ra = vis_ra
    image_dec = vis_dec
    if 'ra_deg' in image_config:
        image_ra = math.radians(image_config['ra_deg'])
    if 'dec_deg' in image_config:
        image_dec = math.radians(image_config['dec_deg'])
    telescope_model = telescope['path']
    lon = math.radians(telescope['lon_deg'])
    lat = math.radians(telescope['lat_deg'])
    alt = telescope['alt_m']
    size = image_config['size']
    fov = image_config['fov_deg']
    x, y, z = load_station_coords(join(telescope_model, 'layout.txt'))
    x, y, z = convert_enu_to_ecef(x, y, z, lon, lat, alt)
    num_antennas = x.shape[0]
    num_baselines = num_antennas * (num_antennas - 1) / 2

    print('- Making image ...')
    print('  - Description =', image_config['description'])
    print('  - Size        =', size)
    print('  - FoV         =', fov)
    print('  - Type         =', image_config['type'])

    # TODO-BM switch on amp = {model|data}
    amp = vis[image_config['type']]
    uu = vis['uu']
    vv = vis['vv']
    ww = vis['ww']

    # fig = pyplot.figure()
    # ax = fig.add_subplot(111)
    # # ax.plot(uu[::10], vv[::10], '.')
    # uv_dist = (uu[::10]**2 + vv[::10]**2)**0.5
    # ax.plot(uv_dist, numpy.abs(amp[::10]), '+')
    # pyplot.show()

    if image_ra != vis_ra or image_dec != vis_dec:
        print('  - Phase rotating visibilities ..')
        l, m, n = convert_ra_dec_to_relative_directions(image_ra, image_dec,
                                                        vis_ra, vis_dec)
        delta_l = -l
        delta_m = -m
        delta_n = 1.0 - n
        phase = 2.0 * math.pi * (uu * delta_l + vv * delta_m + ww * delta_n)
        amp *= numpy.exp(1.0j * phase)
        uu, vv, ww = generate_baseline_uvw(x, y, z, image_ra, image_dec,
                                           num_dumps, num_baselines,
                                           freq_hz, mjd_start,
                                           dump_time_s)

    t0 = time.time()
    im = make(uu, vv, ww, amp, fov, size)

    print('  - Time taken to generate image = %.1f s' % (time.time() - t0))

    return im
