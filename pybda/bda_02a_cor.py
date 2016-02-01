# -*- coding: utf-8 -*-
"""Module to create corruptions."""

import numpy
import os
from os.path import join
import time
import json
import pickle
import shutil
from numpy import array, int32
import utilities


def get_time_info(ms):
    """."""
    tb.open(ms, nomodify=True)
    times = numpy.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [numpy.min(times), numpy.max(times)]
    num_times = len(times)
    length = time_range[1] - time_range[0]
    dt = length / (num_times - 1)
    return num_times, time_range, length, dt


def get_num_antennas(ms):
    """."""
    tb.open(ms + '/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations


def add_extra_data_columns(ms):
    """Add MODEL_DATA and CORRECTED_DATA columns to a Measurement Set."""
    tb.open(ms, nomodify=False)
    # -------------------------------- Add the MODEL_DATA column and fill
    #                                  with values from the 'DATA' column.
    model_desc = {
        'MODEL_DATA': {
            'comment': 'model data',
            'dataManagerGroup': 'ModelTiled',
            'dataManagerType': 'TiledShapeStMan',
            'maxlen': 0,
            'ndim': 2,
            'option': 0,
            'valueType': 'complex'
        }
    }
    model_dminfo = {
        '*6': {
            'COLUMNS': array(['MODEL_DATA'], dtype='|S11'),
            'NAME': 'ModelTiled',
            'SEQNR': 5,
            'SPEC': {
                'ActualMaxCacheSize': 0,
                'DEFAULTTILESHAPE': array([1, 1, 64262], dtype=int32),
                'HYPERCUBES': {
                    '*1': {
                        'BucketSize': 514096,
                        'CellShape': array([1, 1], dtype=int32),
                        'CubeShape': array([1, 1, 1606550], dtype=int32),
                        'ID': {},
                        'TileShape': array([1, 1, 64262], dtype=int32)
                    }
                },
                'IndexSize': 1,
                'MAXIMUMCACHESIZE': 0,
                'SEQNR': 5
            },
            'TYPE': 'TiledShapeStMan'
        }
    }
    tb.addcols(model_desc, model_dminfo)
    tb.putcol('MODEL_DATA', tb.getcol('DATA'))
    # -------------------------------- Add the CORRECTED_DATA column
    #                                  and fill with zeros.
    cor_desc = {
        'CORRECTED_DATA':
        {
            'comment': 'corrected data',
            'dataManagerGroup': 'CorrectedTiled',
            'dataManagerType': 'TiledShapeStMan',
            'maxlen': 0,
            'ndim': 2,
            'option': 0,
            'valueType': 'complex'
        }}
    cor_dminfo = {
        '*7': {
            'COLUMNS': array(['CORRECTED_DATA'], dtype='|S15'),
            'NAME': 'CorrectedTiled',
            'SEQNR': 6,
            'SPEC': {
                'ActualMaxCacheSize': 0,
                'DEFAULTTILESHAPE': array([1, 1, 64262], dtype=int32),
                'HYPERCUBES': {
                    '*1': {
                        'BucketSize': 514096,
                        'CellShape': array([1, 1], dtype=int32),
                        'CubeShape': array([1, 1, 1606550], dtype=int32),
                        'ID': {},
                        'TileShape': array([1, 1, 64262], dtype=int32)
                    }
                },
                'IndexSize': 1,
                'MAXIMUMCACHESIZE': 0,
                'SEQNR': 6
            },
            'TYPE': 'TiledShapeStMan'
        }
    }
    tb.addcols(cor_desc, cor_dminfo)
    tb.putcol('CORRECTED_DATA', numpy.zeros((1, 1, tb.nrows()), dtype='c16'))
    tb.close()


# TODO(BM) Split generation and apply of gains.
def corrupt_data(ms, settings):
    """."""
    num_antennas = get_num_antennas(ms)
    num_times, time_range, length, dt = get_time_info(ms)

    # ----------------------------------
    tau = round(settings['tau_s'] / dt) * dt
    amp_std_t0 = settings['amplitude']['std_t0']
    amp_H = settings['amplitude']['hurst']
    amp_adev_fbm = settings['amplitude']['allan_var']
    amp_sigma_wn = 0.0
    phase_std_t0 = settings['phase']['std_t0']
    phase_H = settings['phase']['hurst']
    phase_adev_fbm = settings['phase']['allan_var']
    phase_sigma_wn = 0.0
    numpy.random.seed(settings['seed'])
    # ----------------------------------
    print '-' * 60
    print 'dt  = %f' % dt
    print 'tau = %f' % tau
    print '-' * 60

    tb.open(ms, nomodify=True)
    Vpq = tb.getcol('DATA')
    tb.close()
    Vpq = Vpq.flatten()
    num_vis = Vpq.shape[0]
    all_gains = {}
    all_gains[0] = numpy.ones((num_times,), dtype='c16')

    inv_gain = numpy.zeros((num_antennas, num_times), dtype='c16')
    inv_conj_gain = numpy.zeros((num_antennas, num_times), dtype='c16')

    # set value for ref. antenna
    inv_gain[0, :] = numpy.ones((num_times,), dtype='c16')
    inv_conj_gain[0, :] = numpy.ones((num_times,), dtype='c16')

    # set gains for rest of antennas.
    for a in range(1, num_antennas):
        g = utilities.eval_complex_gain(num_times, dt, amp_H, amp_adev_fbm,
                                        amp_sigma_wn, phase_H, phase_adev_fbm,
                                        phase_sigma_wn, amp_std_t0,
                                        phase_std_t0, tau)
        all_gains[a] = g
        inv_gain[a, :] = 1. / g
        inv_conj_gain[a, :] = numpy.conj(1. / g)

    out_dir = os.path.dirname(ms)
    pickle.dump(all_gains, open(os.path.join(out_dir, 'gains.pickle'), 'w'))

    # Apply gains to visibilities.
    idx = 0
    for t in range(0, num_times):
        for p in range(0, num_antennas):
            for q in range(p + 1, num_antennas):
                gp = inv_gain[p, t]
                gq = inv_conj_gain[q, t]
                Vpq[idx] = gp * Vpq[idx] * gq
                idx += 1

    tb.open(ms, nomodify=False)
    tb.putcol('DATA', numpy.reshape(Vpq, (1, 1, num_vis)))
    tb.putcol('CORRECTED_DATA', numpy.reshape(Vpq, (1, 1, num_vis)))
    tb.close()


def main(config_file):
    """."""
    tAll = time.time()

    settings = utilities.byteify(json.load(open(config_file)))
    sim_dir = settings['path']
    settings = settings['corrupt']

    # ---------------------------------------------------
    ms_in = join(sim_dir, settings['input_ms'])
    ms = join(sim_dir, settings['output_ms'])
    cal_table = join(sim_dir, settings['gain_table'])
    # ---------------------------------------------------

    if os.path.isdir(ms):
        shutil.rmtree(ms)
    if os.path.isdir(cal_table):
        shutil.rmtree(cal_table)

    t0 = time.time()
    utilities.copytree(ms_in, ms)
    print '+ Coppied simulation data: %.3fs' % (time.time() - t0)

    t0 = time.time()
    add_extra_data_columns(ms)
    print '+ Added MODEL_DATA and CORRECTED_DATA columns: %.3f s' \
          % (time.time() - t0)

    t0 = time.time()
    corrupt_data(ms, settings)
    print '+ Apply of corruptions: %.3f s' % (time.time() - t0)

    print '+ Total time: %.3f s' % (time.time() - tAll)


if __name__ == "__main__":
    os.nice(19)
    main(config_file)

