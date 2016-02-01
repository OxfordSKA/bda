# -*- coding: utf-8 -*-
"""Calibrate with gaincal."""

import os
from os.path import join
import time
from pybda import utilities
import shutil
import json


def run_gaincal(ms, cal_table):
    """
    minblperant => minimum baselines per antenna to solve (default=4)
    minsnr => reject solutions below this snr (default =3.0
    """
    gaincal(vis=ms, caltable=cal_table, field='', spw='', selectdata=False,
            solint='int', refant='0', minblperant=1, minsnr=3.0, gaintype='G',
            calmode='ap', append=False)


def run_applycal(ms, cal_table):
    """."""
    # TODO(BM) replace with own function tested in test_apply_gains.py
    applycal(vis=ms, field='', spw='', selectdata=False,
             gaintable=[cal_table], gainfield=[''], interp=['nearest'],
             calwt=[False], applymode='calonly', flagbackup=False)


if __name__ == "__main__":
    settings = utilities.byteify(json.load(open(config_file)))

    # -------------------------------------------------------------------------
    sim_dir = settings['path']
    calib = settings['calibration']
    overwrite = False
    # -------------------------------------------------------------------------

    # corrupted.ms
    # corrupted_bda.ms
    # corrupted_bda_expanded.ms
    # corrupted_noisy.ms
    # corrupted_noisy_bda.ms
    # corrupted_noisy_bda_expanded.ms
    # (SLOW) corrupted_sub_sampled.ms

    corrupted_root = settings['ms_name']['corrupted']
    calibrated_root = settings['ms_name']['calibrated']
    suffix = settings['ms_modifier']
    ms_names = [
        '',
        '_%s' % (suffix['bda']),
        '_%s_%s' % (suffix['bda'], suffix['expanded']),
        '_%s' % (suffix['noisy']),
        '_%s_%s' % (suffix['noisy'], suffix['bda']),
        '_%s_%s_%s' % (suffix['noisy'], suffix['bda'], suffix['expanded'])
        #'_%s' % (suffix['sub_sampled'])
    ]

    for i, ms in enumerate(ms_names):

        ms_in = join(sim_dir, '%s%s.ms' % (corrupted_root, ms))
        ms_out = join(sim_dir, '%s%s.ms' % (calibrated_root, ms))
        cal_table = join(sim_dir, '%s%s.gains' % (calibrated_root, ms))

        if os.path.isdir(ms_out) or not os.path.isdir(ms_in):
            continue

        print '=' * 80
        print '[%02i/%02i]' % (i + 1, len(ms_names))
        print 'Calibrating:', ms_in
        print '        -->:', ms_out
        print '        -->:', cal_table

        shutil.copytree(ms_in, ms_out)

        t0 = time.time()
        run_gaincal(ms_out, cal_table)
        print '- Gaincal completed in %.3fs' % (time.time() - t0)

        if not os.path.isdir(cal_table):
            print '$' * 80
            print 'ERROR: failed to calibrate!', ms_out
            print '$' * 80
            shutil.rmtree(ms_out)
            continue

        t0 = time.time()
        run_applycal(ms_out, cal_table)
        print '- Applycal completed in %.3fs' % (time.time() - t0)
