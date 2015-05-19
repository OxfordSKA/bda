#!/usr/bin/python

# To ignore undefined CASA classes
#     pylint: disable=undefined-variable

import os
import numpy as np
import time
import shutil


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def load_gains(cal_table):
    tb.open(cal_table, nomodify=True)
    gains = tb.getcol('CPARAM')
    tb.close()
    tb.open(cal_table+'/ANTENNA')
    num_antennas = tb.nrows()
    tb.close()
    num_times = gains.shape[2]/num_antennas
    tb.open(cal_table+'/OBSERVATION')
    time_range = tb.getcol('TIME_RANGE')
    tb.close()
    dt = (time_range[1][0]-time_range[0][0]) / num_times
    return gains, num_antennas, num_times, dt


def load_vis_amp(ms):
    tb.open(ms, nomodify=True)
    data = tb.getcol('DATA')
    tb.close()
    return data


def update_vis_amp(ms, data):
    tb.open(ms, nomodify=False)
    tb.putcol('DATA', data)
    tb.flush()
    tb.close()


def main():
    tAll = time.time()
    cal_table = os.path.join('vis', 'test.cal')
    ms_in = os.path.join('vis', 'test.ms')
    ms = os.path.join('vis', 'test_cor2.ms')

    t0 = time.time()
    if os.path.isdir(ms):
        shutil.rmtree(ms)
    copytree(ms_in, ms)
    print 'Create MS copy = %.3fs' % (time.time() - t0)

    t0 = time.time()
    gains, num_antennas, num_times, dt = load_gains(cal_table)
    data = load_vis_amp(ms)
    print 'Time taken to load data = %.3fs' % (time.time() - t0)

    t0 = time.time()
    gains = gains[0, 0, :]
    g_inv = 1./gains
    g_inv_conj = np.conj(g_inv)
    print 'Time taken to prepare gains = %.3fs' % (time.time() - t0)

    data_shape = data.shape
    data = data.flatten()
    data_out = np.zeros(data.shape, dtype=data.dtype)

    t0 = time.time()
    num_baselines = (num_antennas * (num_antennas-1)) / 2
    i = 0
    for t in range(0, num_times):
        for p in range(0, num_antennas):
            for q in range(p+1, num_antennas):
                ip = t * num_antennas + p
                iq = t * num_antennas + q
                data_out[i] = g_inv[p] * data[i] * g_inv_conj[q]
                #print i, t, p, q, ip, iq
                i += 1
    print 'Time taken to apply gains = %.3fs' % (time.time() - t0)

    t0 = time.time()
    update_vis_amp(ms, np.reshape(data_out, data_shape))
    print 'Time taken to write gains = %.3fs' % (time.time() - t0)

    print 'Total time = %.3fs' % (time.time() - tAll)

if __name__ == "__main__":
    main()
