# -*- coding: utf-8 -*-

import os
import numpy as np
import pickle
import progressbar
import matplotlib.pyplot as plt

def load_ms(ms):
    dtype = [('a1','i4'), ('a2', 'i4'), ('data', 'c16'), ('t', 'f8')]
    tb.open(ms)
    num_rows = tb.nrows()
    values = np.zeros((num_rows,), dtype=dtype)
    values['data'] = np.squeeze(tb.getcol('DATA'))
    values['a1'] = np.squeeze(tb.getcol('ANTENNA1'))
    values['a2'] = np.squeeze(tb.getcol('ANTENNA2'))
    time = tb.getcol('TIME')
    time -= time[0]
    values['t'] = time
    tb.close()
    return values


def load_ms_2(ms):
    dtype = [('a1','i4'), ('a2', 'i4'), ('data', 'c16'), ('t', 'f8'),
             ('model_data', 'c16'), ('corrected_data', 'c16')]
    tb.open(ms)
    num_rows = tb.nrows()
    values = np.zeros((num_rows,), dtype=dtype)
    values['data'] = np.squeeze(tb.getcol('DATA'))
    values['model_data'] = np.squeeze(tb.getcol('MODEL_DATA'))
    values['corrected_data'] = np.squeeze(tb.getcol('CORRECTED_DATA'))
    values['a1'] = np.squeeze(tb.getcol('ANTENNA1'))
    values['a2'] = np.squeeze(tb.getcol('ANTENNA2'))
    time = tb.getcol('TIME')
    time -= time[0]
    values['t'] = time
    tb.close()
    return values


test_model = False
test_corrupted = True
if test_model:
    print 'Model:'
    default_ms = os.path.join('out_default', 'vis', 'model_mstransform.ms')
    new_ms = os.path.join('out_new', 'vis', 'model_bda.ms')

    default = load_ms(default_ms)
    new = load_ms(new_ms)

    print '*' * 60
    print default.shape
    print new.shape
    print '*' * 60

    print 'diffs:'
    print 'data', np.max(default['data'] - new['data'])
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

    print 'diffs:'
    print 'data           :', np.max(np.abs(default['data'] - new['data']))
    print 'model_data     :', np.max(np.abs(default['model_data'] - new['model_data']))
    print 'corrected_data :', np.max(np.abs(default['corrected_data'] - new['corrected_data']))

