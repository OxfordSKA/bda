#!/usr/bin/python
"""Module to create corruptions."""

import numpy as np
import os
import shutil
import time
from numpy import array, int32
from numpy.fft import fft
from numpy.random import randn
import matplotlib.pyplot as plt


def adev(data, dt, tau):
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
                      phase_H, phase_adev_fbm, phase_sigma_wn,
                      amp_std_t0, phase_std_t0, tau):
    amp = eval_gain_amp(n, dt, amp_H, amp_adev_fbm, amp_sigma_wn, tau)
    phase = eval_gain_phase(n, dt, phase_H, phase_adev_fbm, phase_sigma_wn,
                            tau)*np.pi/180.0
    amp_t0 = np.random.randn(1,)*amp_std_t0
    phase_t0 = np.random.randn(1,)*phase_std_t0*(np.pi/180.)
    gain = (amp+amp_t0) * np.exp(1.0j * (phase+phase_t0))
    return gain[0:n]


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


# TODO(BM) split this function into generating and applying gains. The apply
# can then be reused.
def corrupt_data(ms):
    num_antennas = get_num_antennas(ms)
    num_times, time_range, length, dt = get_time_info(ms)

    # ----------------------------------
    tau = np.ceil(1.0/dt)*dt
    amp_std_t0 = 0.005
    amp_H = 0.8
    amp_adev_fbm = 1.0e-4
    amp_sigma_wn = 0.0
    phase_std_t0 = 30.0
    phase_H = 0.8
    phase_adev_fbm = 1.0
    phase_sigma_wn = 0.0

    # ----------------------------------

    # Antenna 0 is the ref antenna so has gain of 1+0j
    # print 'num_times = ', num_times
    # g = eval_complex_gain(num_times, dt, amp_H, amp_adev_fbm, amp_sigma_wn,
    #                       phase_H, phase_adev_fbm, phase_sigma_wn,
    #                       amp_std_t0, phase_std_t0, tau)
    # print g.shape
    # fig = plt.figure()
    # ax = fig.add_subplot(211)
    # ax.plot(np.arange(0, g.shape[0])*dt, np.angle(g)*(180./np.pi), '+:')
    # ax.set_ylabel('Phase [degrees]')
    # ax = fig.add_subplot(212)
    # ax.plot(np.arange(0, g.shape[0])*dt, np.abs(g), '+:')
    # ax.set_ylabel('Amplitude')
    # ax.set_xlabel('Time [seconds]')
    # plt.show()

    # Load visibility amplitudes.
    # dtype=[('a1', 'i4'), ('a2', 'i4'), ('amp', 'c16')]
    # tb.open(ms, nomodify=True)
    # vis = np.zeros((tb.nrows(),), dtype=dtype)
    # vis['a1'] = tb.getcol('ANTENNA1')
    # vis['a2'] = tb.getcol('ANTENNA2')
    # vis['amp'] = tb.getcol('DATA')
    # tb.close()
    # Vpq = vis['amp']
    tb.open(ms, nomodify=True)
    Vpq = tb.getcol('DATA')
    tb.close()
    Vpq = Vpq.flatten()
    num_vis = Vpq.shape[0]

    gain = np.zeros((num_antennas, num_times), dtype='c16')
    conj_gain = np.zeros((num_antennas, num_times), dtype='c16')
    gain[0, :] = np.ones((num_times,),dtype='c16')
    conj_gain[0, :] = np.ones((num_times,),dtype='c16')
    for a in range(1, num_antennas):
        g = eval_complex_gain(num_times, dt, amp_H, amp_adev_fbm, amp_sigma_wn,
                              phase_H, phase_adev_fbm, phase_sigma_wn,
                              amp_std_t0, phase_std_t0, tau)
        g = 1./g
        gain[a, :] = g
        conj_gain[a, :] = np.conj(g)



    # Apply gains to visibilities.
    idx = 0
    for p in range(0, num_antennas):
        gp = gain[p, :]
        for q in range(p+1, num_antennas):
            # Vpq = vis[vis['a1']==p]
            # Vpq = Vpq[Vpq['a2']==q]
            gq = conj_gain[q, :]
            for t in range(0, num_times):
                Vpq[idx] = gp[t]*Vpq[idx]*gq[t]
                idx += 1

    tb.open(ms, nomodify=False)
    tb.putcol('DATA', np.reshape(Vpq, (1, 1, num_vis)))
    tb.close()


def fill_caltable(cal_table, num_stations, num_times, time_range, dt):
    # ----------------------------------
    tau = np.ceil(1.0/dt)*dt
    amp_H = 0.8
    amp_adev_fbm = 1.0e-4
    amp_sigma_wn = 0.0
    phase_H = 0.8
    phase_adev_fbm = 1.0
    phase_sigma_wn = 0.0
    # TODO(BM) generate station allan variance within a distribution so some
    #          stations are better than others wrt their allan variance.
    #          - randomise the gain and phase at the start time for each
    #            antenna. gain = gain0 + gain
    # amp_std_t0 = 0.005
    # phase_std_t0 = 90.0
    # TODO(BM) generate gains wrt ref. antenna
    # ----------------------------------
    tb.open(cal_table, nomodify=False)
    gains0 = eval_complex_gain(num_times, dt,
                               amp_H, amp_adev_fbm, amp_sigma_wn,
                               phase_H, phase_adev_fbm, phase_sigma_wn, tau)
    for t in range(0, num_times):
        row = 0 + t * num_stations
        gain = tb.getcell('CPARAM', row)
        # gain[0] = gains0[t]
        # gain[1] = gains0[t]
        gain[0] = 1.0 + 0.0j
        gain[1] = 1.0 + 0.0j
        tb.putcell('CPARAM', row, gain)
        tb.putcell('TIME', row, time_range[0] + t * dt)

    # TODO(BM) Consider sorting the MS data into baseline time order first?

    tb.open(cal_table, nomodify=False)
    for s in range(1, num_stations):
        # amp_t0 = np.random.randn(1,)*amp_std_t0
        # phase_t0 = np.random.randn(1,)*phase_std_t0
        # gain_t0 = amp_t0 * np.exp(1.0j * phase_t0 * (np.pi/180.0))
        #print amp_t0, phase_t0, np.angle(gain_t0)*(180.0/np.pi)
        gains = eval_complex_gain(num_times, dt,
                                  amp_H, amp_adev_fbm, amp_sigma_wn,
                                  phase_H, phase_adev_fbm, phase_sigma_wn, tau)
        # print np.angle(gains[0])*(180./np.pi)
        #gains += gain_t0
        # print np.angle(gains[0])*(180./np.pi), np.angle(gain_t0)*(180./np.pi)
        # print ''
        #gains = gains/gains0
        for t in range(0, num_times):
            row = s + t * num_stations
            gain = tb.getcell('CPARAM', row)
            gain[0] = gains[t]
            gain[1] = gains[t]
            tb.putcell('CPARAM', row, gain)
            tb.putcell('TIME', row, time_range[0] + t * dt)

    tb.close()



def add_extra_data_columns(ms):
    tb.open(ms, nomodify=False)
    # Add the model data column and fill with values from the 'DATA' column.
    model_desc = {
        'MODEL_DATA': {
            'comment': 'model data',
            'dataManagerGroup': 'ModelTiled',
            'dataManagerType': 'TiledShapeStMan',
            'maxlen': 0,
            'ndim': 2,
            'option': 0,
            'valueType': 'complex'
        }
    }
    model_dminfo = {
        '*6': {
            'COLUMNS': array(['MODEL_DATA'], dtype='|S11'),
            'NAME': 'ModelTiled',
            'SEQNR': 5,
            'SPEC': {
                'ActualMaxCacheSize': 0,
                'DEFAULTTILESHAPE': array([1, 1, 64262], dtype=int32),
                'HYPERCUBES': {
                    '*1': {
                        'BucketSize': 514096,
                        'CellShape': array([1, 1], dtype=int32),
                        'CubeShape': array([1, 1, 1606550], dtype=int32),
                        'ID': {},
                        'TileShape': array([1, 1, 64262], dtype=int32)
                    }
                },
                'IndexSize': 1,
                'MAXIMUMCACHESIZE': 0,
                'SEQNR': 5
            },
            'TYPE': 'TiledShapeStMan'
        }
    }
    tb.addcols(model_desc, model_dminfo)
    tb.putcol('MODEL_DATA', tb.getcol('DATA'))
    # Add the corrected data column and initialise to zero.
    cor_desc = {
        'CORRECTED_DATA': {
            'comment': 'corrected data',
            'dataManagerGroup': 'CorrectedTiled',
            'dataManagerType': 'TiledShapeStMan',
            'maxlen': 0,
            'ndim': 2,
            'option': 0,
            'valueType': 'complex'
            }
        }
    cor_dminfo = {
        '*7': {
            'COLUMNS': array(['CORRECTED_DATA'], dtype='|S15'),
            'NAME': 'CorrectedTiled',
            'SEQNR': 6,
            'SPEC': {
                'ActualMaxCacheSize': 0,
                'DEFAULTTILESHAPE': array([1, 1, 64262], dtype=int32),
                'HYPERCUBES': {
                    '*1': {
                        'BucketSize': 514096,
                        'CellShape': array([1, 1], dtype=int32),
                        'CubeShape': array([1, 1, 1606550], dtype=int32),
                        'ID': {},
                        'TileShape': array([1, 1, 64262], dtype=int32)
                    }
                },
                'IndexSize': 1,
                'MAXIMUMCACHESIZE': 0,
                'SEQNR': 6
            },
            'TYPE': 'TiledShapeStMan'
        }
    }
    tb.addcols(cor_desc, cor_dminfo)
    tb.putcol('CORRECTED_DATA', np.zeros((1, 1, tb.nrows()), dtype='c16'))
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
    ms_in = os.path.join('vis', 'model.ms')
    ms = os.path.join('vis', 'corrupted.ms')
    cal_table = os.path.join('vis', 'corrupted.gains')
    # ---------------------------------------------------

    if os.path.isdir(ms):
        shutil.rmtree(ms)
    if os.path.isdir(cal_table):
        shutil.rmtree(cal_table)

    t0 = time.time()
    copytree(ms_in, ms)
    print '+ Coppied simulation data: %.3fs' % (time.time() - t0)

    t0 = time.time()
    add_extra_data_columns(ms)
    print '+ Added MODEL_DATA and CORRECTED_DATA columns: %.3f s' \
          % (time.time() - t0)

    t0 = time.time()
    corrupt_data(ms)
    print '+ Apply of corruptions: %.3f s' % (time.time() - t0)


    # num_times, time_range, obs_length, dt = get_time_info(ms)
    # num_stations = get_num_antennas(ms)
    # print '-'*80
    # print '+ MS           : %s' % ms
    # print '+ Cal table    : %s' % cal_table
    # print '+ No. times    : %i' % num_times
    # print '+ Time range   : %f %f' % (time_range[0], time_range[1])
    # print '+ Obs. length  : %f' % obs_length
    # print '+ Delta t      : %f' % dt
    # print '+ No. stations : %i' % num_stations
    # print '-'*80
    # print ''
    #
    # t0 = time.time()
    # print '+ Creating calibration table ...'
    # create_empty_caltable(ms, cal_table, num_times)
    # print '+ Done [%.3fs].\n' % (time.time() - t0)
    #
    # t0 = time.time()
    # print '+ Filling calibration table'
    # fill_caltable(cal_table, num_stations, num_times, time_range, dt)
    # print '+ Done [%.3fs].\n' % (time.time() - t0)
    #
    # t0 = time.time()
    # print '+ Applying calibration ...'
    # run_applycal(ms, cal_table)
    # print '+ Done [%.3fs].\n' % (time.time() - t0)
    #
    # t0 = time.time()
    # print '+ Updating data column ...'
    # copy_column(ms, 'CORRECTED_DATA', 'DATA')
    # print '+ Done [%.3fs].\n' % (time.time() - t0)
    #
    # # TODO(BM) add white noise (do this before/ after applying gain errors?!)

    print '+ Total time: %.3f s' % (time.time()-tAll)

if __name__ == "__main__":
    """Run with:
    casa --nologger --nogui -c bda_02_cor.py
    """
    main()
