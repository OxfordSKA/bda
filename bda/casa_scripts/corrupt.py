# -*- coding: utf-8 -*-
"""Module to create corruptions."""

import numpy
import shutil
import os
from os.path import join
import time
import pickle
from bda import utilities
import json


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


def fill_caltable(settings, cal_table, num_stations, num_times, time_range, dt):
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

    out_dir = os.path.dirname(cal_table)
    all_gains = {}
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
        gains = utilities.eval_complex_gain(num_times, dt, amp_H, amp_adev_fbm,
                                            amp_sigma_wn, phase_H,
                                            phase_adev_fbm, phase_sigma_wn,
                                            amp_std_t0, phase_std_t0, tau)
        all_gains[s] = gains
        for t in range(0, num_times):
            row = s + t * num_stations
            gain = tb.getcell('CPARAM', row)
            gain[0] = gains[t]
            gain[1] = gains[t]
            tb.putcell('CPARAM', row, gain)
            tb.putcell('TIME', row, time_range[0] + t * dt)

    gains_pickle = join(out_dir, 'sub_sampled_' + settings['gain_table'] +
                        '.pickle')
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
    """."""
    tAll = time.time()

    settings = utilities.byteify(json.load(open(config_file)))
    sim_dir = settings['path']
    settings = settings['corrupt']

    # ---------------------------------------------------
    ms_in = join(sim_dir, 'sub_sampled_' + settings['input_ms'])
    ms_out = join(sim_dir, 'sub_sampled_' + settings['output_ms'])
    cal_table = join(sim_dir, 'sub_sampled_' + settings['gain_table'])
    # ---------------------------------------------------

    if os.path.isdir(ms_out):
        return

    if os.path.isdir(cal_table):
        shutil.rmtree(cal_table)

    t0 = time.time()
    print '+ Coping simulated MS %s to %s ...' % (ms_in, ms_out)
    utilities.copytree(ms_in, ms_out)
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
    fill_caltable(settings, cal_table, num_stations, num_times, time_range, dt)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Applying corruptions ...'
    run_applycal(ms_out, cal_table)
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    t0 = time.time()
    print '+ Updating data column ...'
    copy_column(ms_out, 'CORRECTED_DATA', 'DATA')
    print '+ Done [%.3fs].\n' % (time.time() - t0)

    # TODO(BM) add white noise (do this before/ after applying gain errors?!)
    print '+ Applying corruptions took %.3f seconds' % (time.time() - tAll)


if __name__ == "__main__":
    os.nice(19)
    main(config_file)
