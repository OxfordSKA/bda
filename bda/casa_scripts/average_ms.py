# -*- coding: utf-8 -*-

import numpy
import os
from os.path import join
import shutil
import json
from bda import utilities


def _get_num_antennas(ms):
    tb.open(ms + '/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations


def _average_column(data, num_blocks, num_times_per_block, num_baselines):
    data = data.reshape(num_blocks, num_times_per_block, num_baselines)
    ave_data = numpy.mean(data, axis=1)
    ave_data = ave_data.reshape(num_blocks * num_baselines)
    return ave_data


def _average_ms(ms_ref, ms_in, ms_out, num_baselines, num_times_average,
               overwrite=True):
    if not overwrite and os.path.isdir(ms_out):
        return
    # --- Create output MS by making a copy of the reference MS.
    if os.path.exists(ms_out):
        shutil.rmtree(ms_out)
    print 'Averaging MS:', ms_in
    shutil.copytree(ms_ref, ms_out)
    # --- Open the input MS and average.
    tb.open(ms_in, nomodify=True)
    num_rows = tb.nrows()
    num_times = num_rows / num_baselines
    num_blocks = num_times / num_times_average
    # TODO-BM: average CORRECTED_DATA AND MODEL_DATA columns as well..

    col_data = tb.getcol('DATA')
    # col_uvw = tb.getcol('UVW')
    # col_ant1 = tb.getcol('ANTENNA1')
    # col_ant2 = tb.getcol('ANTENNA2')
    # col_time = tb.getcol('TIME')

    # ave_uu = average_column(col_uvw[0, :], num_blocks, num_times_average,
    #                         num_baselines)
    # ave_vv = average_column(col_uvw[1, :], num_blocks, num_times_average,
    #                         num_baselines)
    # ave_ww = average_column(col_uvw[2, :], num_blocks, num_times_average,
    #                         num_baselines)
    # ave_t = average_column(col_time, num_blocks, num_times_average,
    #                        num_baselines)
    # Assert that the MS has 1 channel and is stokes-I only.
    assert col_data.shape[0] == 1
    assert col_data.shape[1] == 1
    assert col_data.shape[2] == num_rows
    ave_data = _average_column(numpy.squeeze(col_data), num_blocks,
                               num_times_average, num_baselines)
    if 'MODEL_DATA' in tb.colnames():
        print 'MODEL_DATA found in ', ms_in
        col_model_data = tb.getcol('MODEL_DATA')
        ave_model_data = _average_column(numpy.squeeze(col_model_data),
                                         num_blocks, num_times_average,
                                         num_baselines)
    if 'CORRECTED_DATA' in tb.colnames():
        print 'CORRECTED_DATA found in ', ms_in
        col_corrected_data = tb.getcol('CORRECTED_DATA')
        ave_corrected_data = _average_column(numpy.squeeze(col_corrected_data),
                                             num_blocks, num_times_average,
                                             num_baselines)
    tb.close()
    # --- Open the output MS and write averaged values.
    tb.open(ms_out, nomodify=False)
    col_data = tb.getcol('DATA')
    tb.putcol('DATA', numpy.reshape(ave_data, col_data.shape))
    if not 'MODEL_DATA' in tb.colnames():
        if not 'MODEL_DATA' in tb.colnames():
            print '=>> Adding MODEL_DATA column to', ms_out
            tb.close()
            clearcal(vis=ms_out, addmodel=True)
            tb.open(ms_out, nomodify=False)
        tb.putcol('MODEL_DATA', numpy.reshape(ave_model_data,
                                              col_data.shape))
    if not 'CORRECTED_DATA' in tb.colnames():
        if not 'CORRECTED_DATA' in tb.colnames():
            print '=>> Adding CORRECTED_DATA column to', ms_out
            tb.close()
            clearcal(vis=ms_out, addmodel=True)
            tb.open(ms_out, nomodify=False)
        tb.putcol('CORRECTED_DATA', numpy.reshape(ave_corrected_data,
                                                  col_data.shape))
    tb.close()


if __name__ == "__main__":
    """Copy the ref ms and populate it with averaged input ms."""
    settings = utilities.byteify(json.load(open(config_file)))
    sim_dir = settings['path']
    ms_ref = join(sim_dir, 'ref_' + settings['sim']['output_ms'])
    num_antennas = _get_num_antennas(ms_ref)
    num_baselines = num_antennas * (num_antennas - 1) / 2
    num_times_average = settings['sim']['observation']['over_sample']
    overwrite = False

    ms_in = join(sim_dir, 'sub_sampled_' + settings['sim']['output_ms'])
    ms_out = join(sim_dir, settings['sim']['output_ms'])
    _average_ms(ms_ref, ms_in, ms_out, num_baselines, num_times_average,
                overwrite=overwrite)

    ms_in = join(sim_dir, 'sub_sampled_' + settings['corrupt']['output_ms'])
    ms_out = join(sim_dir, settings['corrupt']['output_ms'])
    _average_ms(ms_ref, ms_in, ms_out, num_baselines, num_times_average,
                overwrite=overwrite)
