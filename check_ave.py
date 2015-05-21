#!/usr/bin/python

import os
import numpy as np


dtype=[('t',  'f8'), ('uu', 'f8'), ('vv', 'f8'), ('ww', 'f8'),
        ('a1', 'i4'), ('a2', 'i4'), ('v', 'c16')]
#dtype='f8,f8,f8,f8,i4,i4,c16'


ms = os.path.join('vis', 'corrupted.ms')
tb.open(ms)
orig_data = tb.getcol('DATA')
orig_uvw = tb.getcol('UVW')
orig_time = tb.getcol('TIME')
orig_time_delta = orig_time-orig_time[0]
orig_ant1 = tb.getcol('ANTENNA1')
orig_ant2 = tb.getcol('ANTENNA2')
tb.close()
orig = np.zeros((orig_time.shape[0],), dtype=dtype)
orig['t'] = orig_time_delta
orig['uu'] = orig_uvw[0, :]
orig['vv'] = orig_uvw[1, :]
orig['ww'] = orig_uvw[2, :]
orig['a1'] = orig_ant1
orig['a2'] = orig_ant2
orig['v'] = orig_data[0, 0, :]

p = 0
q = 1
orig_pq = orig[orig['a1']==p]
orig_pq = orig_pq[orig_pq['a2']==q]
print '+ No. vis for baseline_pq in original data %i' % orig_pq.shape[0]
orig_intervals = np.round(orig_pq['t']/1.6)


ms = os.path.join('vis', 'corrupted_ave.ms')
tb.open(ms)
ave_data = tb.getcol('DATA')
ave_uvw = tb.getcol('UVW')
ave_time = tb.getcol('TIME')
ave_time_delta = ave_time-ave_time[0]
ave_ant1 = tb.getcol('ANTENNA1')
ave_ant2 = tb.getcol('ANTENNA2')
tb.close()
ave = np.zeros((ave_time.shape[0],), dtype=dtype)

ave['t'] = ave_time_delta
ave['uu'] = ave_uvw[0, :]
ave['vv'] = ave_uvw[1, :]
ave['ww'] = ave_uvw[2, :]
ave['a1'] = ave_ant1
ave['a2'] = ave_ant2
ave['v'] = ave_data[0, 0, :]
# Obtain baseline 01 visibilities
ave_pq = ave[ave['a1']==p]
ave_pq = ave_pq[ave_pq['a2']==q]
print '+ No. vis for baseline_pq in averaged data %i' % ave_pq.shape[0]
ave_intervals = ave_pq['t']/1.6

# FIXME(BM) mmm not sure why this doenst work...
start=0
for i in range(0, 20):
    diff_pq = ave_pq['v'][0]-(np.mean(orig_pq['v'][start:start+i+1]))
    #diff_pq = ave_pq['v'][0]-orig_pq['v'][start+i]
    print '%02i %+7.3f %+7.3f [%+7.3f] [%+8.3f]' % \
        (i, np.real(diff_pq), np.imag(diff_pq), np.abs(diff_pq),
         np.angle(diff_pq)*(180./np.pi))
print ''

for i in range(0, 20):
    diff_pq_0 = ave_pq['t'][0]/1.6-orig_pq['t'][start+i]/1.6
    diff_pq_1 = ave_pq['t'][1]/1.6-orig_pq['t'][start+i]/1.6
    print '%02i %+7.3f %+7.3f' % (i, diff_pq_0, diff_pq_1)
    if i == 7 or i == 19: print ''


# TODO(BM) try with different ms transform settings.
# ie. perhaps with timebin to dt*2 for testing.

print ''
print orig.shape[0], ave.shape[0], '%.2f:1' % (float(orig.shape[0])/ave.shape[0])
