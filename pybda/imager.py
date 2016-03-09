# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import
import time
import numpy
from os.path import join
from oskar.imager import Imager

def run(config, vis, amp_name, image_name=None):
    images = []
    border_trim = 0.4 # Fraction of 1
    for image_id in range(len(config['imaging']['images'])):
        sim_dir = config['path']
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
        num_vis = len(vis[amp_name])

        print('- Making image...')
        print('  - Description            :', image_config['description'])
        print('  - Amplitude array name   :', amp_name)
        print('  - Number of visibilities :', num_vis)
        print('  - Size                   :', size)
        print('  - FoV                    : %.1f deg' % fov)

        # Set up the imager.
        im = numpy.zeros([size, size], dtype='f8')
        img = Imager("double")
        img.set_fov(fov)
        img.set_size(size)
        if image_name != None:
            img.set_output_root(join(sim_dir, image_name + '_' + str(image_id)))
        img.set_vis_frequency(freq_hz, 0.0, 1)
        img.set_vis_time(mjd_start + 0.5 * (num_dumps * dump_time_s), 
            dump_time_s, 1)
        img.set_vis_phase_centre(vis_ra, vis_dec)
        if image_ra != vis_ra or image_dec != vis_dec:
            img.set_direction(image_ra, image_dec)
        block_size = num_vis / 100
        num_blocks = (num_vis + block_size - 1) / block_size
        t0 = time.time()
        for i in range(num_blocks):
            start = i * block_size
            end = start + block_size
            if end > num_vis:
                end = num_vis
            img.update(end-start, vis['uu'][start:end], vis['vv'][start:end],
                vis['ww'][start:end], vis[amp_name][start:end],
                vis['weight'][start:end])
        img.finalise(im)
        pix_start = border_trim * size
        pix_end = size - pix_start
        images.append(im[pix_start:pix_end, pix_start:pix_end])
        print('  - Visibilities imaged in %.1f s' % (time.time() - t0))

    return images
