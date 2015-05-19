#!/usr/bin/python
"""Module to create corruptions."""

import numpy as np
import os
import shutil
import time

def adev(data, dt, tau):
    """Evaluate the Allan deviation of a time series.
    Args:
        data (np.array): Time series data.
        dt (float): Time series data increment, in seconds.
        tau (float): Averaging period, in seconds.
    Returns:
        The allan deviation of the series and error on the deviation.
    """
    rate = 1./dt
    m = int(rate * tau)
    # Truncate to an even multiple of this tau value
    freq = data[0:len(data)-int(np.remainder(len(data), m))]
    f = np.reshape(freq, (m, -1), order='F')
    fa = np.mean(f, 0)
    fd = np.diff(fa)
    M = len(fa)
    M = M-1
    sm = np.sqrt(0.5/(M)*(np.sum(fd**2)))
    sme = sm / np.sqrt(M)
    return sm, sme, M


def fbm(n, H, seed=None):
    """Generate Fractional brownian motion noise.

    http://www.maths.uq.edu.au/~kroese/ps/MCSpatial.pdf

    Args:
        n (int): Length of the time series.
        H (float): Hurst parameter.

    Kwargs:
        seed (int): Random geerator number seed.

    Returns:
        Time series array.
    """
    from numpy.fft import fft
    from numpy.random import randn
    if seed:
        np.random.seed(seed)  # Seed the random number generator
    r = np.empty((n+1,))*np.NaN
    r[0] = 1
    for k in range(1, n+1):
        a = 2.0*H
        r[k] = 0.5*((k+1)**a - 2*k**a + (k-1)**a)
    r = np.append(r, r[-2:0:-1])  # first row of circulant matrix
    lambda_ = np.real(fft(r))/(2*n)  # Eigenvalues
    W = fft(np.sqrt(lambda_)*(randn(2*n,)+1j*randn(2*n,)))
    W = n**(-H)*np.cumsum(np.real(W[0:n+1]))  # Rescale
    return W


def eval_gain_amp(length, tstep, hurst, adev_fbm, sigma_wn, tau):
    g_fbm = fbm(length, hurst)
    g_fbm *= adev_fbm/adev(g_fbm, tstep, tau)[0]
    g_wn = np.random.randn(length+1,)*sigma_wn
    return (g_fbm + g_wn) + 1.0

def eval_gain_phase(n, tstep, H, adev_fbm, sigma_wn, tau):
    p_fbm = fbm(n, H)
    p_fbm *= adev_fbm/adev(p_fbm, tstep, tau)[0]
    p_wn = np.random.randn(n+1,)*sigma_wn
    phase = p_fbm + p_wn  # Combine FBM and WN signals
    return phase


def eval_complex_gain(n, dt, amp_H, amp_adev_fbm, amp_sigma_wn,
                      phase_H, phase_adev_fbm, phase_sigma_wn, tau):
    amp = eval_gain_amp(n, dt, amp_H, amp_adev_fbm, amp_sigma_wn, tau)
    phase = eval_gain_phase(n, dt, phase_H, phase_adev_fbm, phase_sigma_wn,
                            tau)*np.pi/180.0
    gain = amp * np.exp(1j * phase)
    return gain


def get_time_info(ms):
    tb.open(ms, nomodify=True)
    times = np.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [np.min(times), np.max(times)]
    num_times = len(times)
    length = time_range[1]-time_range[0]
    dt = length/(num_times-1)
    return num_times, time_range, length, dt


def get_num_antennas(ms):
    tb.open(ms+'/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations


def create_empty_caltable(ms, cal_table, num_times):
    if os.path.isdir(cal_table):
        shutil.rmtree(cal_table)
    # Gencal generates a table with one row for each antenna
    gencal(vis=ms, caltable=cal_table, caltype='amp', parameter=[1.0])
    # Copy all rows for each remaining time interval (num_times-1)
    tb.open(cal_table, nomodify=False)
    num_rows = tb.nrows()
    for i in range(0, num_times-1):
        tb.copyrows(outtable=cal_table, startrowin=0, startrowout=-1,
                    nrow=num_rows)
    tb.close()


def fill_caltable(cal_table, num_stations, num_times, time_range, dt):
    # ----------------------------------
    tau = np.ceil(1.0/dt)*dt
    amp_H = 0.8
    amp_adev_fbm = 2.0e-3
    amp_sigma_wn = 0.0
    phase_H = 0.8
    phase_adev_fbm = 0.15
    phase_sigma_wn = 0.0
    # TODO(BM) generate station allan variance within a distribution so some
    #          stations are better than others wrt their allan variance.
    # ----------------------------------
    tb.open(cal_table, nomodify=False)
    for s in range(0, num_stations):
        gains = eval_complex_gain(num_times, dt,
                                  amp_H, amp_adev_fbm, amp_sigma_wn,
                                  phase_H, phase_adev_fbm, phase_sigma_wn,
                                  tau)
        for t in range(0, num_times):
            row = s + t * num_stations
            gain = tb.getcell('CPARAM', row)
            gain[0] = gains[t]
            gain[1] = gains[t]
            tb.putcell('CPARAM', row, gain)
            tb.putcell('TIME', row, time_range[0] + t * dt)
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
    if not 'CORRECTED_DATA' in colnames or not 'MODEL_DATA' in colnames:
        clearcal(vis=ms, addmodel=True)


def copy_column(ms, src_col, dst_col):
    tb.open(ms, nomodify=False)
    num_rows = tb.nrows()
    rows_read = 0
    while rows_read != num_rows:
        chunk_size = 10000
        if num_rows - rows_read < chunk_size:
            chunk_size = num_rows - rows_read
        t = tb.getcol(src_col, startrow=rows_read, nrow=chunk_size)
        tb.putcol(dst_col, value=t, startrow=rows_read, nrow=chunk_size)
        rows_read += chunk_size
    tb.close()

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def main():
    tAll = time.time()

    # ---------------------------------------------------
    ms_in = os.path.join('vis', 'test.ms')
    ms = os.path.join('vis', 'test_cor.ms')
    cal_table = os.path.join('vis', 'test.cal')
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
    print '-'*80
    print '+ MS           : %s' % ms
    print '+ Cal table    : %s' % cal_table
    print '+ No. times    : %i' % num_times
    print '+ Time range   : %f %f' % (time_range[0], time_range[1])
    print '+ Obs. length  : %f' % obs_length
    print '+ Delta t      : %f' % dt
    print '+ No. stations : %i' % num_stations
    print '-'*80
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

    print '+ Applying corruptions took %.3f seconds' % (time.time()-tAll)

if __name__ == "__main__":
    """Run with:
    casa --nologger --nogui -c bda_02_cor.py
    """
    main()
