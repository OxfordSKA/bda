# -*- coding: utf-8 -*-
"""Module to create corruptions."""

import numpy as np
import sys
import os
import time
import pickle
sys.path.append(os.path.join(os.getcwd(), 'scripts'))
from bda_00_util import *


def get_time_info(ms):
    """."""
    tb.open(ms, nomodify=True)
    times = np.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [np.min(times), np.max(times)]
    num_times = len(times)
    length = time_range[1] - time_range[0]
    dt = length / (num_times - 1)
    return num_times, time_range, length, dt


def get_num_antennas(ms):
    """."""
    tb.open(ms + '/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations


def create_empty_caltable(ms, cal_table, num_times):
    if os.path.isdir(cal_table):
        shutil.rmtree(cal_table)
    # gencal() generates a table with one row for each antenna.
    gencal(vis=ms, caltable=cal_table, caltype='amp', parameter=[1.0])
    # Copy all rows for each remaining time interval (num_times-1).
    tb.open(cal_table, nomodify=False)
    num_rows = tb.nrows()
    for i in range(0, num_times - 1):
        tb.copyrows(outtable=cal_table, startrowin=0, startrowout=-1,
                    nrow=num_rows)
    tb.close()


def fill_caltable(cal_table, num_stations, num_times, time_range, dt):
    """Fill a CASA calibration table."""
    # ----------------------------------
    tau = np.ceil(1.0 / dt)*dt
    amp_std_t0 = 0.005
    amp_H = 0.8
    amp_adev_fbm = 1.0e-4
    amp_sigma_wn = 0.0
    phase_std_t0 = 30.0
    phase_H = 0.8
    phase_adev_fbm = 1.0
    phase_sigma_wn = 0.0
    np.random.seed(101)
    # ----------------------------------
    out_dir = os.path.dirname(cal_table)
    all_gains = {}
    all_gains[0] = np.ones((num_times,), dtype='c16')
    tb.open(cal_table, nomodify=False)
    for t in range(0, num_times):
        row = 0 + t * num_stations
        gain = tb.getcell('CPARAM', row)
        gain[0] = 1.0 + 0.0j
        gain[1] = 1.0 + 0.0j
        tb.putcell('CPARAM', row, gain)
        tb.putcell('TIME', row, time_range[0] + t * dt)

    tb.open(cal_table, nomodify=False)
    for s in range(1, num_stations):
        gains = eval_complex_gain(num_times, dt,
                                  amp_H, amp_adev_fbm, amp_sigma_wn,
                                  phase_H, phase_adev_fbm, phase_sigma_wn,
                                  amp_std_t0, phase_std_t0, tau)
        all_gains[s] = gains
        for t in range(0, num_times):
            row = s + t * num_stations
            gain = tb.getcell('CPARAM', row)
            gain[0] = gains[t]
            gain[1] = gains[t]
            tb.putcell('CPARAM', row, gain)
            tb.putcell('TIME', row, time_range[0] + t * dt)

    pickle.dump(all_gains, open(os.path.join(out_dir, 'gains.pickle'), 'w'))

    tb.close()


def run_applycal(ms, cal_table):
    # Performs Vis_cal = 1./G_p * V_pq * conj(1./G_q)
    applycal(vis=ms, field='', spw='', selectdata=False, gaintable=[cal_table],
             gainfield=[''], interp=['nearest'], calwt=[False],
             applymode='calonly', flagbackup=False)


def scratch_columns_create(ms):
    tb.open(ms)
    colnames = tb.colnames()
    tb.close()
    if 'CORRECTED_DATA' not in colnames or 'MODEL_DATA' not in colnames:
        clearcal(vis=ms, addmodel=True)


def copy_column(ms, src_col, dst_col):
    tb.open(ms, nomodify=False)
    num_rows = tb.nrows()
    rows_read = 0
    while rows_read != num_rows:
        chunk_size = 20000
        if num_rows - rows_read < chunk_size:
            chunk_size = num_rows - rows_read
        t = tb.getcol(src_col, startrow=rows_read, nrow=chunk_size)
        tb.putcol(dst_col, value=t, startrow=rows_read, nrow=chunk_size)
        rows_read += chunk_size
    tb.close()


def main(sim_dir):
    tAll = time.time()

    # ---------------------------------------------------
    ms_in = os.path.join(sim_dir, 'vis', 'model.ms')
    ms = os.path.join(sim_dir, 'vis', 'corrupted.ms')
    cal_table = os.path.join(sim_dir, 'vis', 'corrupted.gains')
    # ---------------------------------------------------

    if os.path.isdir(ms):
        shutil.rmtree(ms)
    if os.path.isdir(cal_table):
        shutil.rmtree(cal_table)

    t0 = time.time()
    print '+ Coping simulated MS %s to %s ...' % (ms_in, ms)
    copytree(ms_in, ms)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Creating scratch (CORRECTED_DATA & MODEL_DATA) columns in MS ...'
    scratch_columns_create(ms)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Copy DATA to MODEL_DATA ...'
    copy_column(ms, 'DATA', 'MODEL_DATA')
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    num_times, time_range, obs_length, dt = get_time_info(ms)
    num_stations = get_num_antennas(ms)
    print '-' * 80
    print '+ MS           : %s' % ms
    print '+ Cal table    : %s' % cal_table
    print '+ No. times    : %i' % num_times
    print '+ Time range   : %f %f' % (time_range[0], time_range[1])
    print '+ Obs. length  : %f' % obs_length
    print '+ Delta t      : %f' % dt
    print '+ No. stations : %i' % num_stations
    print '-' * 80
    print ''

    t0 = time.time()
    print '+ Creating calibration table ...'
    create_empty_caltable(ms, cal_table, num_times)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Filling calibration table'
    fill_caltable(cal_table, num_stations, num_times, time_range, dt)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Applying calibration ...'
    run_applycal(ms, cal_table)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Updating data column ...'
    copy_column(ms, 'CORRECTED_DATA', 'DATA')
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    # TODO(BM) add white noise (do this before/ after applying gain errors?!)
    print '+ Applying corruptions took %.3f seconds' % (time.time() - tAll)


if __name__ == "__main__":
    if len(sys.argv) - 1 < 1:
        print 'Usage:'
        print ('  $ casa --nologger --nogui -c scripts/bda_02_cor.py '
               '<simulation dir>')
        sys.exit(1)

    sim_dir = sys.argv[-1]
    if not os.path.isdir(sim_dir):
        print 'ERROR: simulation directory not found!'
        sys.exit(1)

    print '-' * 60
    print 'Simulation directory:', sim_dir
    print '-' * 60

    os.nice(19)
    main(sim_dir)
