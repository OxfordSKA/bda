#!/usr/bin/python

import os
import numpy as np


def inv_sinc(arg):
    """
    Newton-Raphson method for calculating arcsinc(x), from Obit.
    """
    import numpy as np
    x1 = 0.001
    for i in range(0, 1000):
        x0 = x1
        a = x0 * np.pi
        x1 = x0-((np.sin(a)/a)-arg) / ((a*np.cos(a) - np.pi*np.sin(a)) / (a**2))
        if (np.fabs(x1 - x0) < 1.0e-6): break
    return x1


dtype=[('t',  'f8'), ('uu', 'f8'), ('vv', 'f8'), ('ww', 'f8'),
        ('a1', 'i4'), ('a2', 'i4'), ('int', 'f8'), ('v', 'c16')]

# msroot = 'model'
msroot = 'corrupted'
ms = os.path.join('test_sim_01', 'vis', '%s.ms' % msroot)
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


ms = os.path.join('test_sim_01', 'vis', '%s_ave.ms' % msroot)
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

# Obtain baseline 01 visibilities
p = 10
q = 50
orig_pq = orig[orig['a1']==p]
orig_pq = orig_pq[orig_pq['a2']==q]
ave_pq = ave[ave['a1']==p]
ave_pq = ave_pq[ave_pq['a2']==q]

print ''
print 'No. vis orig = %i' % orig.shape[0]
print 'No. vis ave  = %i' % ave.shape[0]
print 'Compression  = %.2f:1' % (float(orig.shape[0])/ave.shape[0])
print ''

print '+ Looking at baseline pq (%i,%i)' % (p,q)
print '+ No. vis for baseline_pq in original data %i' % orig_pq.shape[0]
print '+ No. vis for baseline_pq in averaged data %i' % ave_pq.shape[0]
print ''

# Looking at the times long baseline pq ...
# dt is likely not 1.6 due to precision drop in storing as MJD seconds.
dt0 = orig_pq['t'][1]-orig_pq['t'][0]
print 'dt0 = %.15f' % dt0
dt_orig = orig_pq['t'][0]-orig_pq['t'][0]
diff = ave_pq['t'][0]-(orig_pq['t'][0]-dt0/2.0)
print 't0 dt_orig %+7.3f d_diff: ave-orig: %+7.3f' % (dt_orig/dt0, diff/dt0)
print ''

# Try to evalute the average along baseline pq ...
# Seems to be averaging on in blocks of 10 as mean(orig[0:10]) ~ ave[0]
print '# diff = ave_pq[0] - mean(orig_pq[0:i+1])'
print '# i  real(diff) imag(diff) [abs(diff)] [angle(diff)]'
for i in range(7, 12):
    diff_pq = ave_pq['v'][0]-(np.mean(orig_pq['v'][0:i+1]))
    print '%02i %+7.3e %+7.3e [%+7.3e] [%+8.3e]' % \
        (i, np.real(diff_pq), np.imag(diff_pq), np.abs(diff_pq),
         np.angle(diff_pq)*(180./np.pi))
print ''

# uv coordinates (basic test to check if uvw coordinates are averaged or
# reevaluated for new time...)
# TODO(BM) check this over all baselines and times.
# TODO(BM) wipe out antenna table to prevent recalculation and see if BDA
# still works.
print 'ave_uvw[0] =', ave_pq['uu'][0], ave_pq['vv'][0], ave_pq['ww'][0]
uu_ = np.mean(orig_pq['uu'][0:10])
vv_ = np.mean(orig_pq['vv'][0:10])
ww_ = np.mean(orig_pq['ww'][0:10])
print 'mean(orig_uvw[0:10]) =', uu_, vv_, ww_
print 'diff: ave_uvw[0]-mean(orig_uvw[0:10]) =', ave_pq['uu'][0]-uu_,\
                                                ave_pq['vv'][0]-vv_,\
                                                ave_pq['ww'][0]-ww_

# TODO(BM) work out averaging intervals from
#   timebin(dt_max) and
#   maxuvwdistance(delta_uv)
#
max_fact = 1.01
freq = 700.0e6
fov_radius = 0.9
wavelength = 299792458.0 / freq
duvw_max = inv_sinc(1.0 / max_fact) / (fov_radius * (np.pi/180.0))
duvw_max *= wavelength
dt_max = 10.0*dt0

print ''
print '* wavelength = %f m' % wavelength
print '* delta_uv   = %f m' % duvw_max
print '* dt_max     = %f seconds' % dt_max
print ''

num_times = len(orig_pq['t'])
i0 = 0
for i in range(0, num_times):
    duu = orig_pq['uu'][i]-orig_pq['uu'][i0]
    dvv = orig_pq['vv'][i]-orig_pq['vv'][i0]
    dww = orig_pq['ww'][i]-orig_pq['ww'][i0]
    duvw = (duu**2+dvv**2+dww**2)*0.5
    dt = (orig_pq['t'][i]-orig_pq['t'][i0])
    luvw = (orig_pq['uu'][i]**2+orig_pq['vv'][i]**2+orig_pq['ww'][i]**2)**0.5
    print '%2i -  %f  - %f - %f' %  (i, duvw, luvw, dt)
    if (duvw > duvw_max or dt > dt_max):
        print '*** %i' % (i-i0+1)
        i0 = i+1
print '*** %i' % (i-i0+1)
