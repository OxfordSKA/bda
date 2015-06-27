# -*- coding: utf-8 -*-
"""Calibrate with gaincal."""

import os
import time
import sys
sys.path.append(os.path.join(os.getcwd(), 'scripts'))
from bda_00_util import *


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


def main(sim_dir):
    """."""
    non_bda = True
    bda = True

    if non_bda:
        # -------------------------------------------------------------------------
        ms_in = os.path.join(sim_dir, 'vis', 'corrupted.ms')
        ms = os.path.join(sim_dir, 'vis', 'calibrated.ms')
        cal_table = os.path.splitext(ms)[0] + '.gains'
        # -------------------------------------------------------------------------

        if os.path.isdir(ms):
            print 'Removing existing MS : %s' % ms
            shutil.rmtree(ms)
        os.system('cp -r %s %s' % (ms_in, ms))
        # copytree(ms_in, ms)

        t0 = time.time()
        run_gaincal(ms, cal_table)
        print '*' * 80
        print '+ Gaincal completed in %.3fs' % (time.time() - t0)
        print '*' * 80
        print '\n' * 3

        t0 = time.time()
        run_applycal(ms, cal_table)
        print '*' * 80
        print '+ Applycal completed in %.3fs' % (time.time() - t0)
        print '*' * 80
        print '\n' * 3

    if bda:
        # -------------------------------------------------------------------------
        ms_in = os.path.join(sim_dir, 'vis', 'corrupted_mstransform.ms')
        ms = os.path.join(sim_dir, 'vis', 'calibrated_mstransform.ms')
        if not os.path.isdir(ms_in):
            ms_in = os.path.join(sim_dir, 'vis', 'corrupted_bda.ms')
        if not os.path.isdir(ms):
            ms = os.path.join(sim_dir, 'vis', 'calibrated_bda.ms')
        cal_table = os.path.splitext(ms)[0] + '.gains'
        # -------------------------------------------------------------------------

        if os.path.isdir(ms):
            print 'Removing existing MS : %s' % ms
            shutil.rmtree(ms)

        # copytree(ms_in, ms)
        os.system('cp -r %s %s' % (ms_in, ms))

        t0 = time.time()
        run_gaincal(ms, cal_table)
        print '*' * 80
        print '+ Gaincal completed in %.3fs' % (time.time() - t0)
        print '*' * 80
        print '\n' * 3

        t0 = time.time()
        run_applycal(ms, cal_table)
        print '*' * 80
        print '+ Applycal completed in %.3fs' % (time.time() - t0)
        print '*' * 80


if __name__ == "__main__":
    if len(sys.argv) - 1 < 1:
        print 'Usage:'
        print ('  $ casa --nologger --nogui -c scripts/bda_04_cal.py '
               '<simulation dir>')
        sys.exit(1)

    sim_dir = sys.argv[-1]
    if not os.path.isdir(sim_dir):
        print 'ERROR: simulation directory not found!'
        sys.exit(1)

    main(sim_dir)
