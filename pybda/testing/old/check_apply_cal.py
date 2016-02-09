#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as pp
import time

def load_gains(cal_table):
    dtype=[('time',  'f8'), ('a1', 'i4'), ('a2', 'i4'), ('cparam', 'c16', 2),
           ('paramerr', 'f8', 2), ('flag', 'i4',2), ('snr', 'f8', 2)]
    tb.open(cal_table, nomodify=True)
    gains = np.zeros((tb.nrows(),), dtype=dtype)
    gains['time'] = tb.getcol('TIME')
    gains['a1'] = tb.getcol('ANTENNA1')
    gains['a2'] = tb.getcol('ANTENNA2')
    gains_cparam = tb.getcol('CPARAM')
    gains['cparam'][:, 0] = gains_cparam[0, 0, :]
    gains['cparam'][:, 1] = gains_cparam[1, 0, :]
    gains_paramerr = tb.getcol('PARAMERR')
    gains['paramerr'][:, 0] = gains_paramerr[0, 0, :]
    gains['paramerr'][:, 1] = gains_paramerr[1, 0, :]
    gains_flag = tb.getcol('FLAG')
    gains['flag'][:, 0] = gains_flag[0, 0, :]
    gains['flag'][:, 1] = gains_flag[1, 0, :]
    gains_snr = tb.getcol('SNR')
    gains['snr'][:, 0] = gains_snr[0, 0, :]
    gains['snr'][:, 1] = gains_snr[1, 0, :]
    tb.close()

    tb.open(cal_table+'/ANTENNA')
    num_antennas = tb.nrows()
    tb.close()

    return gains, num_antennas


def load_ms(ms, data_column='DATA'):
    dtype=[('time',  'f8'), ('uu', 'f8'), ('vv', 'f8'), ('ww', 'f8'),
           ('a1', 'i4'), ('a2', 'i4'), ('interval', 'f8'), ('vis', 'c16')]
    tb.open(ms)
    vis = np.zeros((tb.nrows(),), dtype=dtype)
    vis['time'] = tb.getcol('TIME')
    uvw_ = tb.getcol('UVW')
    vis['uu'] = uvw_[0, :]
    vis['vv'] = uvw_[1, :]
    vis['ww'] = uvw_[2, :]
    vis['a1'] = tb.getcol('ANTENNA1')
    vis['a2'] = tb.getcol('ANTENNA2')
    vis['interval'] = tb.getcol('INTERVAL')
    data_ = tb.getcol(data_column)
    vis['vis'] = data_[0, 0, :]
    tb.close()
    return vis


load_data = True
if load_data:
    # ----------------------------------------------------------------------
    cal_table = os.path.join('vis', 'calibrated_ave.gains')
    cor_ms = os.path.join('vis', 'corrupted_ave.ms')
    cal_ms = os.path.join('vis', 'calibrated_ave.ms')
    ave_ms = os.path.join('vis', 'model_ave.ms')
    # ----------------------------------------------------------------------

    t0 = time.time()
    cal_gains, num_antennas = load_gains(cal_table)
    cor_vis = load_ms(cor_ms)
    cal_vis = load_ms(cal_ms, 'CORRECTED_DATA')
    ave_vis = load_ms(ave_ms)
    print 'Time taken to load data = %.3fs' % (time.time() - t0)
 
 
p = 30
q = 225
gains_p = cal_gains[cal_gains['a1']==p]
gains_p = gains_p[gains_p['a2']==0]
gains_q = cal_gains[cal_gains['a1']==q]
gains_q = gains_q[gains_q['a2']==0]

cor_vis_pq = cor_vis[cor_vis['a1']==p]
cor_vis_pq = cor_vis_pq[cor_vis_pq['a2']==q]
cal_vis_pq = cal_vis[cal_vis['a1']==p]
cal_vis_pq = cal_vis_pq[cal_vis_pq['a2']==q]
ave_vis_pq = ave_vis[ave_vis['a1']==p]
ave_vis_pq = ave_vis_pq[ave_vis_pq['a2']==q]


# print gains_p['cparam'].shape
# print gains_q['cparam'].shape
# print cal_vis_pq['vis'].shape
# print cor_vis_pq['vis'].shape
print ''

plot_gains = True
if plot_gains:
    fig = pp.figure(1)
    ax = fig.add_subplot(111)
    x = gains_q['time']-gains_q['time'][0]
    yp = np.abs(gains_p['cparam'][:,0])
    yq = np.abs(gains_q['cparam'][:,0])
    ax.plot(x, yp, 'ro--')
    ax.plot(x, 1./np.conj(yq), 'bo--')
    #y = yp*1./np.conj(yq)
    #ax.plot(x, y, 'go--')
    pp.show()

print 'Gains:'
diff = []
for i in range(0, gains_p.shape[0]):
    gp = gains_p['cparam'][i,0]
    gq = gains_q['cparam'][i,0]
    diff.append((1.0/gp * cor_vis_pq['vis'][0] * 1.0/np.conj(gq))- cal_vis_pq['vis'][0])
    print diff[-1]

fig = pp.figure(2)
ax = fig.add_subplot(111)
x = gains_q['time'] # -gains_q['time'][0]
y = np.abs(diff)
ax.plot(x, y, 'ro--')
ax.plot([cor_vis_pq['time'][0], cor_vis_pq['time'][0]], [0, np.max(np.abs(diff))], lw=2.0)
ax.plot([cor_vis_pq['time'][1], cor_vis_pq['time'][1]], [0, np.max(np.abs(diff))], lw=2.0)
pp.show()
    

print ''
print 'Vis:'
for i in range(0, cal_vis_pq['vis'].shape[0]):
    m = ave_vis_pq['vis'][i]
    v = cor_vis_pq['vis'][i]
    c = cal_vis_pq['vis'][i]
    print i, v, c



