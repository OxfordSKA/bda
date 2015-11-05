# -*- coding: utf-8 -*-
"""Module to create corruptions."""

import numpy
import shutil
import os
from os.path import join
import time
import pickle


from bda.utilities import eval_complex_gain, byteify
import json


def _smooth(x, window_len=11, window='hanning'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")
    s = numpy.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='valid')
    return y


def get_time_info(ms):
    """."""
    tb.open(ms, nomodify=True)
    times = numpy.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [numpy.min(times), numpy.max(times)]
    num_times = len(times)
    length = time_range[1] - time_range[0]
    dt = length / (num_times - 1)
    return num_times, time_range, length + dt, dt


def get_num_antennas(ms):
    """."""
    tb.open(ms + '/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations


def create_empty_caltable(ms, cal_table, num_times):
    """."""
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


def fill_caltable(settings, cal_table, num_stations, num_times, time_range, dt,
                  over_sample):
    """Fill a CASA calibration table."""
    # ----------------------------------
    tau = round(settings['tau_s'] / dt) * dt
    amp_std_t0 = settings['amplitude']['std_t0']
    amp_H = settings['amplitude']['hurst']
    amp_adev_fbm = settings['amplitude']['allan_var']
    amp_sigma_wn = 0.0
    phase_std_t0 = settings['phase']['std_t0']
    phase_H = settings['phase']['hurst']
    phase_adev_fbm = settings['phase']['allan_var']
    phase_sigma_wn = 0.0
    numpy.random.seed(settings['seed'])
    # ----------------------------------
    print '-' * 60
    print 'dt  = %f' % dt
    print 'tau = %f' % tau
    print 'obs. length = %f' % ((time_range[1] - time_range[0]) + dt)
    print '-' * 60

    smooth_gains = settings['smooth']
    window_length = settings['smooth_length']
    if smooth_gains:
        print 'INFO: Smoothing gains with window length = %i' % window_length
    # TODO-BM 1. set all gains to 1+0j
    # TODO-BM 2. vary all gains no ref antenna and include smoothing (perhaps on the window length of the oversample?)
    all_gains = dict()
    all_gains[0] = numpy.ones((num_times,), dtype='c16')
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
        gains = eval_complex_gain(num_times, dt, amp_H, amp_adev_fbm,
                                  amp_sigma_wn, phase_H,
                                  phase_adev_fbm, phase_sigma_wn,
                                  amp_std_t0, phase_std_t0, tau)
        if smooth_gains:
            gains = _smooth(gains, window_len=window_length)
            gains = gains[0:num_times]
        all_gains[s] = gains
        # all_gains[s] = numpy.ones((num_times,), dtype='c16')

        for t in range(0, num_times):
            row = s + t * num_stations
            gain = tb.getcell('CPARAM', row)
            gain[0] = all_gains[s][t]
            gain[1] = all_gains[s][t]
            tb.putcell('CPARAM', row, gain)
            tb.putcell('TIME', row, time_range[0] + t * dt)

    gains_pickle = '%s.pickle' % cal_table
    pickle.dump(all_gains, open(gains_pickle, 'w'))
    tb.close()


def run_applycal(ms, cal_table):
    """Apply a gain table.

    Performs Vis_cal = 1./G_p * V_pq * conj(1./G_q)
    """
    applycal(vis=ms, field='', spw='', selectdata=False, gaintable=[cal_table],
             gainfield=[''], interp=['nearest'], calwt=[False],
             applymode='calonly', flagbackup=False)


def scratch_columns_create(ms):
    """Adds scratch columns (CORRECTED_DATA and MODEL_DATA) to the ms."""
    tb.open(ms)
    colnames = tb.colnames()
    tb.close()
    if 'CORRECTED_DATA' not in colnames or 'MODEL_DATA' not in colnames:
        clearcal(vis=ms, addmodel=True)


def copy_column(ms, src_col, dst_col):
    """Copy data from one column to another in a Measurement Set."""
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


def main(config_file):
    tAll = time.time()
    settings = byteify(json.load(open(config_file)))
    sim_dir = settings['path']

    # -------------------------------------------------------------------------
    ms_in = join(sim_dir,
                 '%s_%s.ms' % (settings['ms_name']['model'],
                               settings['ms_modifier']['sub_sampled']))
    ms_out = join(sim_dir,
                  '%s_%s.ms' % (settings['ms_name']['corrupted'],
                                settings['ms_modifier']['sub_sampled']))
    cal_table = join(sim_dir,
                     '%s_%s.gains' % (settings['ms_name']['corrupted'],
                                      settings['ms_modifier']['sub_sampled']))
    # -------------------------------------------------------------------------

    if os.path.isdir(ms_out):
        return

    if os.path.isdir(cal_table):
        shutil.rmtree(cal_table)

    t0 = time.time()
    print '+ Coping simulated MS %s to %s ...' % (ms_in, ms_out)
    shutil.copytree(ms_in, ms_out)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Creating scratch (CORRECTED_DATA & MODEL_DATA) columns in MS ...'
    scratch_columns_create(ms_out)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Copy DATA to MODEL_DATA ...'
    copy_column(ms_out, 'DATA', 'MODEL_DATA')
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    num_times, time_range, obs_length, dt = get_time_info(ms_out)
    num_stations = get_num_antennas(ms_out)
    print '-' * 80
    print '+ MS           : %s' % ms_out
    print '+ Cal table    : %s' % cal_table
    print '+ No. times    : %i' % num_times
    print '+ Time range   : %f %f' % (time_range[0], time_range[1])
    print '+ Obs. length  : %f' % obs_length
    print '+ Delta t      : %f' % dt
    print '+ No. stations : %i' % num_stations
    print '-' * 80
    print ''

    t0 = time.time()
    print '+ Creating calibration table for corruptions ...'
    create_empty_caltable(ms_out, cal_table, num_times)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Filling calibration table with corruptions ...'
    fill_caltable(settings['corrupt'], cal_table, num_stations, num_times,
                  time_range, dt, settings['sim']['observation']['over_sample'])
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Applying corruptions ...'
    run_applycal(ms_out, cal_table)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Updating data column ...'
    copy_column(ms_out, 'CORRECTED_DATA', 'DATA')
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    print '+ Applying corruptions took %.3f seconds' % (time.time() - tAll)


if __name__ == "__main__":
    os.nice(19)
    main(config_file)
