# -*- coding: utf-8 -*-
"""Calibration module for the BDA pipeline"""

from pybda.stefcal import stefcal1
from oskar._bda_utils import vis_list_to_matrix, apply_gains
import numpy
import time


def run(vis, verbose=False):
    print('- Calibrating data...')
    num_baselines = vis['num_baselines']
    num_antennas = vis['num_antennas']
    vis['corrected'] = numpy.zeros_like(vis['data'])
    t0 = time.time()
    for i in range(vis['num_times']):
        i0 = i * num_baselines
        i1 = i0 + num_baselines
        data = vis_list_to_matrix(vis['data'][i0:i1], num_antennas)
        model = vis_list_to_matrix(vis['model'][i0:i1], num_antennas)
        g, nit, dg = stefcal1(data, model)
        gains = 1.0 / g
        vis['corrected'][i0:i1] = apply_gains(vis['data'][i0:i1], gains)
        if verbose:
            print(i, nit, dg)
    print('  - Calibration took %.2f s' % (time.time() - t0))

