#!/usr/bin/python -u
"""Check BDA between CASA mstransform() and alternatives."""

import os
import numpy as np
from progressbar import ProgressBar, Percentage, Bar, ETA


def get_delta_uv_max(ms, fov_radius, max_fact):
    """Evaluate maximum allowed baseline coordinate momvement.

    Evalaute the maximum uvw distance allowed to achieve a the specified
    amplitude drop for the specified radius.
    """
    tb.open(ms + '/SPECTRAL_WINDOW', nomodify=True)
    freqs = tb.getcell('CHAN_FREQ', 0)
    tb.close()
    if len(freqs) != 1:
        print 'ERROR: can only use single channel data.'
        return
    freq = freqs[0]
    wavelength = 299792458.0 / freq
    delta_uv = inv_sinc(1.0 / max_fact) / (fov_radius * (np.pi / 180.))
    delta_uv *= wavelength
    return delta_uv


def inv_sinc(arg):
    """Newton-Raphson method for calculating arcsinc(x), from Obit."""
    x1 = 0.001
    for i in range(0, 1000):
        x0 = x1
        a = x0 * np.pi
        x1 = x0 - ((np.sin(a) / a) - arg) / \
            ((a * np.cos(a) - np.pi * np.sin(a)) / (a**2))
        if (np.fabs(x1 - x0) < 1.0e-6):
            break
    return x1


def get_num_antennas(ms):
    """Return number of antennas in a MS."""
    tb.open(ms + '/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations


def load_ms(ms):
    """Load a MS."""
    print '+ Loading data from : %s' % ms
    tb.open(ms)
    dtype = [('time', 'f8'),
             ('ctime', 'f8'),
             ('uvw', 'f8'), ('uu', 'f8'), ('vv', 'f8'),
             ('ww', 'f8'), ('a1', 'i4'), ('a2', 'i4'), ('interval', 'f8'),
             ('exposure', 'f8'), ('data', 'c16')]
    if 'MODEL_DATA' in tb.colnames():
        dtype.append(('model', 'c16'))
    if 'CORRECTED_DATA' in tb.colnames():
        dtype.append(('corrected', 'c16'))
    vis = np.zeros((tb.nrows(),), dtype=dtype)
    vis['time'] = tb.getcol('TIME')
    vis['ctime'] = tb.getcol('TIME_CENTROID')
    uvw_ = tb.getcol('UVW')
    vis['uu'] = uvw_[0, :]
    vis['vv'] = uvw_[1, :]
    vis['ww'] = uvw_[2, :]
    vis['a1'] = tb.getcol('ANTENNA1')
    vis['a2'] = tb.getcol('ANTENNA2')
    vis['interval'] = tb.getcol('INTERVAL')
    vis['exposure'] = tb.getcol('EXPOSURE')
    vis['data'] = tb.getcol('DATA')
    if 'MODEL_DATA' in tb.colnames():
        vis['model'] = tb.getcol('MODEL_DATA')
    if 'CORRECTED_DATA' in tb.colnames():
        vis['corrected'] = tb.getcol('CORRECTED_DATA')
    tb.close()
    return vis


# -----------------------------------------------------------------------------
ms = os.path.join('vis', 'model.ms')
ms_bda_1 = os.path.join('vis', 'model_mstransform_ave.ms')
ms_bda_2 = os.path.join('vis', 'model_bda.ms')
# -----------------------------------------------------------------------------

num_antennas = get_num_antennas(ms)

# Load measurement sets.
if os.path.isdir(ms):
    vis = load_ms(ms)
if os.path.isdir(ms_bda_1):
    vis_bda_1 = load_ms(ms_bda_1)
if os.path.isdir(ms_bda_2):
    vis_bda_2 = load_ms(ms_bda_2)


# Loop over all baselines to find averaging mismatches.
check_fails = True
max_diff = 0
if check_fails:
    print 'No. rows in %s : %i' % (ms_bda_1, vis_bda_1.shape[0])
    print 'No. rows in %s : %i' % (ms_bda_2, vis_bda_2.shape[0])
    print 'Checking MS for averaging differences...'
    widgets = [Percentage(), ' ', Bar(), ' ~', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=len(vis))
    pbar.start()
    i = 0
    for p in range(0, num_antennas):
        for q in range(p + 1, num_antennas):
            pbar.update(i + 1)
            vis_bda_1_pq = vis_bda_1[vis_bda_1['a1'] == p]
            vis_bda_1_pq = vis_bda_1_pq[vis_bda_1_pq['a2'] == q]
            vis_bda_2_pq = vis_bda_2[vis_bda_2['a1'] == p]
            vis_bda_2_pq = vis_bda_2_pq[vis_bda_2_pq['a2'] == q]
            diff = len(vis_bda_1_pq) - len(vis_bda_2_pq)
            max_diff = max(max_diff, diff)
            i += 1
            if diff:
                print '%-3i %-3i : %i %s' % \
                      (p, q, diff, ('***' if diff else ''))
    pbar.finish()
    print 'max diff : %i' % max_diff

# Check averaging condition for a single baseline
else:
    p = 0
    # q = 196
    # q = 168
    q = 199
    vis_pq = vis[vis['a1'] == p]
    vis_pq = vis_pq[vis_pq['a2'] == q]
    vis_bda_1_pq = vis_bda_1[vis_bda_1['a1'] == p]
    vis_bda_1_pq = vis_bda_1_pq[vis_bda_1_pq['a2'] == q]
    vis_bda_2_pq = vis_bda_2[vis_bda_2['a1'] == p]
    vis_bda_2_pq = vis_bda_2_pq[vis_bda_2_pq['a2'] == q]

    print ''
    print 'Baseline : %i, %i' % (p, q)
    print ''
    print '%-30s : %-8i %i' % (ms, vis.shape[0], vis_pq.shape[0])
    print '%-30s : %-8i %i' % (ms_bda_1, vis_bda_1.shape[0],
                               vis_bda_1_pq.shape[0])
    print '%-30s : %-8i %i' % (ms_bda_2, vis_bda_2.shape[0],
                               vis_bda_2_pq.shape[0])
    print ''

    max_fact = 1.01
    fov_deg = 0.9
    dt_max = 5.0 * 1.6
    duv_max = get_delta_uv_max(ms, fov_deg, max_fact)
    print 'duv_max = %.3f' % duv_max
    print 'dt_max  = %.3f' % dt_max

    print 'baseline [%i,%i] no. times = %i' % (p, q, len(vis_pq))
    print 'mstransform no. times = %i' % (len(vis_bda_1_pq))
    print 'custom bda no. times  = %i' % (len(vis_bda_2_pq))
    print ''
    print 't, sum_uu, ave_uu, uu, duvw, duv'
    print '-' * 40
    bi0 = 0
    bt0 = vis_pq['time'][0]
    buu0 = vis_pq['uu'][0]
    bvv0 = vis_pq['vv'][0]
    bww0 = vis_pq['ww'][0]
    sum_uu = 0.0
    count = 0
    for i in range(0, 10):
        count += 1
        dt = vis_pq['time'][i] - bt0 + 1.6
        duu = vis_pq['uu'][i] - buu0
        dvv = vis_pq['vv'][i] - bvv0
        dww = vis_pq['ww'][i] - bww0

        duvw = (duu**2 + dvv**2 + dww**2)**0.5
        # print duvw, count,
        duvw += duvw / float(count - 1)
        # print duvw

        ave_uu = np.mean(vis_pq['uu'][bi0:i + 1])  # array range is exclusive
        sum_uu += vis_pq['uu'][i]

        print '%-2i uu:%12.5f (%9.2f) duvw:%.3f dt:%.3f' % \
              (i, ave_uu, sum_uu, duvw, dt),
        print 'time' if dt >= dt_max else 'uvdist' if duvw >= duv_max else ''
        # If we should move to a new average do so and print /  write current
        # average
        if (dt >= dt_max or duvw >= duv_max or i == 9):
            print '** %i (%i:%i, %i) %f %f' % \
                (i, bi0, i + 1, count,
                 np.mean(vis_pq['uu'][bi0:i + 1]),
                 sum_uu / count)
            print ''
            if i < 9:
                bi0 = i + 1
                bt0 = vis_pq['time'][i + 1]
                buu0 = vis_pq['uu'][i + 1]
                bvv0 = vis_pq['vv'][i + 1]
                bww0 = vis_pq['ww'][i + 1]
                sum_uu = 0.0
                count = 0

    print ''

    print 'Averaging from mstransform:'
    for i in range(0, vis_bda_1_pq.shape[0]):
        print '  %i %f %f : %f%+fj' % (i,
                                       vis_bda_1_pq['uu'][i],
                                       vis_bda_1_pq['time'][i],
                                       np.real(vis_bda_1_pq['data'][i]),
                                       np.imag(vis_bda_1_pq['data'][i]))

    print ''
    print 'Averaging from custom BDA:'
    for i in range(0, vis_bda_2_pq.shape[0]):
        print '  %i %f %f : %f%+fj' % (i,
                                       vis_bda_2_pq['uu'][i],
                                       vis_bda_2_pq['time'][i],
                                       np.real(vis_bda_2_pq['data'][i]),
                                       np.imag(vis_bda_2_pq['data'][i]))
