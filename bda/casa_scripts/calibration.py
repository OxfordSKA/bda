# -*- coding: utf-8 -*-
"""Calibrate with gaincal."""

import os
from os.path import join
import time
from bda import utilities
import shutil
import json


def run_gaincal(ms, cal_table):
    """."""
    gaincal(vis=ms, caltable=cal_table, field='', spw='', selectdata=False,
            solint='int', refant='0', minblperant=1, minsnr=3, gaintype='G',
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

    ms_files = [f for f in os.listdir(os.path.abspath(sim_dir))
                if f.endswith('.ms') and os.path.isdir(join(sim_dir, f))]

    for ms in ms_files:
        if settings['corrupt']['output_ms'] in ms:
            ms_in = join(sim_dir, ms)
            prefix = ''
            for key in settings['ms_prefix']:
                if settings['ms_prefix'][key] in ms:
                    prefix = settings['ms_prefix'][key]
            # if prefix == settings['ms_prefix']['sub_sampled']:
            #     continue
            cal_table = join(sim_dir,
                             prefix + settings['calibration']['output_ms'] +
                             '.gains')
            ms_out = join(sim_dir, prefix +
                          settings['calibration']['output_ms'])
            if not overwrite and os.path.isdir(ms_out):
                continue
            if os.path.isdir(ms_out):
                shutil.rmtree(ms_out)
            if os.path.isdir(cal_table):
                shutil.rmtree(cal_table)
            shutil.copytree(ms_in, ms_out)
            print 'Calibrating:', ms_in
            print '  cal table:', cal_table
            t0 = time.time()
            run_gaincal(ms_out, cal_table)
            print '*' * 80
            print '+ Gaincal completed in %.3fs' % (time.time() - t0)
            print '*' * 80
            t0 = time.time()
            run_applycal(ms_out, cal_table)
            print '*' * 80
            print '+ Applycal completed in %.3fs' % (time.time() - t0)
            print '*' * 80
