# -*- coding: utf-8 -*-

import os
import numpy as np
import pickle
import progressbar
import matplotlib.pyplot as plt


def load_ms(ms):
    dtype = [('uu', 'f8'), ('vv', 'f8'), ('ww', 'f8'),
             ('time', 'f8'), ('time_c', 'f8'),
             ('a1', 'i4'), ('a2', 'i4'), ('data', 'c16'), ('t', 'f8')]
    tb.open(ms)
    num_rows = tb.nrows()
    values = np.zeros((num_rows,), dtype=dtype)
    values['data'] = np.squeeze(tb.getcol('DATA'))
    values['a1'] = np.squeeze(tb.getcol('ANTENNA1'))
    values['a2'] = np.squeeze(tb.getcol('ANTENNA2'))
    uvw = tb.getcol('UVW')
    values['uu'] = uvw[0, :]
    values['vv'] = uvw[1, :]
    values['ww'] = uvw[2, :]
    time = tb.getcol('TIME')
    values['time'] = time
    time -= time[0]
    values['t'] = time
    time_c = tb.getcol('TIME_CENTROID')
    values['time_c'] = time_c
    tb.close()
    return values


def load_ms_2(ms):
    dtype = [('uu', 'f8'), ('vv', 'f8'), ('ww', 'f8'),
             ('time', 'f8'), ('time_c', 'f8'),
             ('a1', 'i4'), ('a2', 'i4'), ('data', 'c16'), ('t', 'f8'),
             ('model_data', 'c16'), ('corrected_data', 'c16')]
    tb.open(ms)
    num_rows = tb.nrows()
    values = np.zeros((num_rows,), dtype=dtype)
    values['data'] = np.squeeze(tb.getcol('DATA'))
    values['model_data'] = np.squeeze(tb.getcol('MODEL_DATA'))
    values['corrected_data'] = np.squeeze(tb.getcol('CORRECTED_DATA'))
    values['a1'] = np.squeeze(tb.getcol('ANTENNA1'))
    values['a2'] = np.squeeze(tb.getcol('ANTENNA2'))
    uvw = tb.getcol('UVW')
    print uvw.shape
    #  values['uu']
    time = tb.getcol('TIME')
    values['time'] = time
    time -= time[0]
    values['t'] = time
    time_c = tb.getcol('TIME_CENTROID')
    values['time_c'] = time_c
    tb.close()
    return values

test_model = True
test_model_ave = True
test_corrupted = True

if test_model:
    print 'Model:'
    default_ms = os.path.join('out_default', 'vis', 'model.ms')
    new_ms = os.path.join('out_new', 'vis', 'model.ms')

    default = load_ms(default_ms)
    new = load_ms(new_ms)

    print '*' * 60
    print default.shape
    print new.shape
    print '*' * 60

    print 'Diffs:'
    print '  DATA  :', np.max(np.abs(default['data'] - new['data']))
    print '  TIME  :', np.max(np.abs(default['time'] - new['time']))
    print '  TIMEC :', np.max(np.abs(default['time_c'] - new['time_c']))
    print '  UU    :', np.max(np.abs(default['uu'] - new['uu']))
    print '  VV    :', np.max(np.abs(default['vv'] - new['vv']))
    print '  WW    :', np.max(np.abs(default['ww'] - new['ww']))
    print '  ANT1  :', np.max(np.abs(default['a1'] - new['a1']))
    print '  ANT2  :', np.max(np.abs(default['a2'] - new['a2']))
    print ''


if test_model_ave:
    print 'Model ave:'
    default_ms = os.path.join('out_default', 'vis', 'model_mstransform.ms')
    new_ms = os.path.join('out_new', 'vis', 'model_bda.ms')

    default = load_ms(default_ms)
    new = load_ms(new_ms)

    print '*' * 60
    print default.shape
    print new.shape
    print '*' * 60

    # diff_time = np.abs(default['time'] - new['time'])
    # plt.plot(diff_time[1250:1300], '+')
    # plt.show()
    # diff_data = np.abs(default['data'] - new['data'])
    # print len(diff_time[diff_time > 1.0])
    # print default['time'][1000]
    # print new['time'][1000]
    # plt.plot(diff_data[0:10000], '+')
    # plt.show()

    print 'Diffs:'
    print '  DATA  :', np.max(np.abs(default['data'] - new['data']))
    print '  TIME  :', np.max(np.abs(default['time'] - new['time']))
    print '  TIMEC :', np.max(np.abs(default['time_c'] - new['time_c']))
    print '  T     :', np.max(np.abs(default['t'] - new['t']))
    print '  UU    :', np.max(np.abs(default['uu'] - new['uu']))
    print '  VV    :', np.max(np.abs(default['vv'] - new['vv']))
    print '  WW    :', np.max(np.abs(default['ww'] - new['ww']))
    print '  ANT1  :', np.max(np.abs(default['a1'] - new['a1']))
    print '  ANT2  :', np.max(np.abs(default['a2'] - new['a2']))
    print ''

if test_corrupted:
    print 'Corrupted:'
    default_ms = os.path.join('out_default', 'vis', 'corrupted_mstransform.ms')
    new_ms = os.path.join('out_new', 'vis', 'corrupted_bda.ms')

    default = load_ms_2(default_ms)
    new = load_ms_2(new_ms)

    print '*' * 60
    print default.shape
    print new.shape
    print '*' * 60

    print 'Diffs:'
    print '  DATA  :', np.max(np.abs(default['data'] - new['data']))
    print '  TIME  :', np.max(np.abs(default['time'] - new['time']))
    print '  TIMEC :', np.max(np.abs(default['time_c'] - new['time_c']))
    print '  T     :', np.max(np.abs(default['t'] - new['t']))
    print '  UU    :', np.max(np.abs(default['uu'] - new['uu']))
    print '  VV    :', np.max(np.abs(default['vv'] - new['vv']))
    print '  WW    :', np.max(np.abs(default['ww'] - new['ww']))
    print '  ANT1  :', np.max(np.abs(default['a1'] - new['a1']))
    print '  ANT2  :', np.max(np.abs(default['a2'] - new['a2']))
    print ''

