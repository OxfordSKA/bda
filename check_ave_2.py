#!/usr/bin/python

import os
import numpy as np


dtype=[('t',  'f8'), ('uu', 'f8'), ('vv', 'f8'), ('ww', 'f8'),
        ('a1', 'i4'), ('a2', 'i4'), ('int', 'f8'), ('v', 'c16')]

msroot = 'model'
ms = os.path.join('vis', '%s.ms' % msroot)
tb.open(ms)
orig_data = tb.getcol('DATA')
orig_uvw = tb.getcol('UVW')
orig_time = tb.getcol('TIME')
orig_ant1 = tb.getcol('ANTENNA1')
orig_ant2 = tb.getcol('ANTENNA2')
orig_int = tb.getcol('INTERVAL')
tb.close()
orig = np.zeros((orig_time.shape[0],), dtype=dtype)
orig['t'] = orig_time
orig['uu'] = orig_uvw[0, :]
orig['vv'] = orig_uvw[1, :]
orig['ww'] = orig_uvw[2, :]
orig['a1'] = orig_ant1
orig['a2'] = orig_ant2
orig['int'] = orig_int
orig['v'] = orig_data[0, 0, :]

orig_01 = orig[orig['a1']==0]
orig_01 = orig_01[orig_01['a2']==1]


ms = os.path.join('vis', '%s_ave.ms' % msroot)
tb.open(ms)
ave_data = tb.getcol('DATA')
ave_uvw = tb.getcol('UVW')
ave_time = tb.getcol('TIME')
ave_ant1 = tb.getcol('ANTENNA1')
ave_ant2 = tb.getcol('ANTENNA2')
ave_int = tb.getcol('INTERVAL')
tb.close()
ave = np.zeros((ave_time.shape[0],), dtype=dtype)
ave['t'] = ave_time
ave['uu'] = ave_uvw[0, :]
ave['vv'] = ave_uvw[1, :]
ave['ww'] = ave_uvw[2, :]
ave['a1'] = ave_ant1
ave['a2'] = ave_ant2
ave['int'] = ave_int
ave['v'] = ave_data[0, 0, :]

ave_01 = ave[ave['a1']==0]
ave_01 = ave_01[ave_01['a2']==1]


# Obtain baseline pq visibilities
max_num_times = 0
for p in range(0, 256):
    for q in range(p+1, 256):
        ave_pq = ave[ave['a1']==p]
        ave_pq = ave_pq[ave_pq['a2']==q]
        max_num_times = max(max_num_times, ave_pq['t'].shape[0])
        print p, q, ave_pq['t'].shape, max_num_times

        
        


