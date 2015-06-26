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


default_ms = os.path.join('out_default', 'vis', 'model_mstransform.ms')
new_ms = os.path.join('out_new', 'vis', 'model_bda.ms')

# default_ms = os.path.join('out_default', 'vis', 'corrupted_mstransform.ms')
# new_ms = os.path.join('out_new', 'vis', 'corrupted_bda.ms')

default = load_ms(default_ms)
new = load_ms(new_ms)

print '*' * 60
print default.shape
print new.shape
print '*' * 60
