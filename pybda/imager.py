# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import
import time
import numpy
from oskar.imager import Imager


def run_imager(config, vis, amp_name, image_name=None):
    images = []
    for image_id in range(len(config['imaging']['images'])):
        image_config = config['imaging']['images'][image_id]
        obs_config = config['sim']['observation']
        mjd_start = obs_config['start_time_mjd']
        num_dumps = obs_config['num_times']
        dump_time_s = obs_config['dump_time_s']
        freq_hz = obs_config['freq_hz']
        vis_ra = obs_config['ra_deg']
        vis_dec = obs_config['dec_deg']
        image_ra = vis_ra
        image_dec = vis_dec
        if 'ra_deg' in image_config:
            image_ra = image_config['ra_deg']
        if 'dec_deg' in image_config:
            image_dec = image_config['dec_deg']
        size = image_config['size']
        fov = image_config['fov_deg']

        amp = vis[amp_name]
        uu = vis['uu']
        vv = vis['vv']
        ww = vis['ww']
        weight = vis['weight']
        num_vis = len(amp)

        print('- Making image...')
        print('  - Description :', image_config['description'])
        print('  - Size        :', size)
        print('  - FoV         : %.1f deg' % fov)
        im = numpy.zeros([size, size], dtype='f8')
        img = Imager("double")
        img.set_fov(fov)
        img.set_size(size)
        #img.set_output_root("output_image_file_name_root") # Enable for FITS file.
        img.set_vis_frequency(freq_hz, 0.0, 1)
        img.set_vis_time(mjd_start + 0.5 * (num_dumps * dump_time_s), 
            dump_time_s, 1)
        img.set_vis_phase_centre(vis_ra, vis_dec)
        if image_ra != vis_ra or image_dec != vis_dec:
            img.set_direction(image_ra, image_dec)
        t0 = time.time()
        img.update(num_vis, uu, vv, ww, amp, weight)
        img.finalise(im)
        images.append(im)
        print('  - Visibilities imaged in %.1f s' % (time.time() - t0))

    return images
