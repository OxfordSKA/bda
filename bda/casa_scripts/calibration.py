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
    # -------------------------------------------------------------------------

    ms_files = [f for f in os.listdir(os.path.abspath(sim_dir))
                if f.endswith('.ms') and os.path.isdir(join(sim_dir, f))]

    for ms in ms_files:
        if settings['corrupt']['output_ms'] in ms:
            print ms

    # for p in zip(settings['input_ms'], settings['output_ms'],
    #              settings['gain_table']):
    #     ms_in = join(sim_dir, p[0])
    #     ms = join(sim_dir, p[1])
    #     cal_table = join(sim_dir, p[2])
    #     if os.path.isdir(ms_in) and not os.path.isdir(ms):
    #         print '-- Calibrating %s' % ms
    #         if os.path.isdir(ms):
    #             print '  * Removing existing MS : %s' % ms
    #             shutil.rmtree(ms)
    #         utilities.copytree(ms_in, ms)
    #         # os.system('cp -r %s %s' % (ms_in, ms))
    #
    #         t0 = time.time()
    #         run_gaincal(ms, cal_table)
    #         print '*' * 80
    #         print '+ Gaincal completed in %.3fs' % (time.time() - t0)
    #         print '*' * 80
    #         print '\n' * 3
    #
    #         t0 = time.time()
    #         run_applycal(ms, cal_table)
    #         print '*' * 80
    #         print '+ Applycal completed in %.3fs' % (time.time() - t0)
    #         print '*' * 80
    #         print '\n' * 3
