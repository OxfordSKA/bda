# -*- coding: utf-8 -*-

from bda import utilities
import os
from os.path import join
import shutil
import json
import numpy


def _get_num_antennas(ms):
    tb.open(ms + '/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations


def run():
    settings = utilities.byteify(json.load(open(config_file)))
    sim_dir = settings['path']

    # We are going to expand the corrupted BDA MS so construct paths for this
    ms_ref = join(sim_dir, settings['ms_prefix']['reference'] +
                  settings['sim']['output_ms'])
    ms_in = join(sim_dir, settings['ms_prefix']['bda'] +
                 settings['corrupt']['output_ms'])
    ms_out = join(sim_dir, settings['ms_prefix']['expanded'] +
                  settings['corrupt']['output_ms'])

    assert os.path.isdir(ms_ref)
    assert os.path.isdir(ms_in)

    # # FIXME-BM remove this delete after testing
    # shutil.rmtree(ms_out)

    if os.path.isdir(ms_out):
        return

    # Copy reference
    shutil.copytree(ms_ref, ms_out)

    # Add model and corrected data columns to the measurement set.
    clearcal(vis=ms_out, addmodel=True)
    # Get output array(s) to write into.
    tb.open(ms_out, nomodify=True)
    out_col_data = tb.getcol('DATA')
    out_col_model_data = tb.getcol('MODEL_DATA')
    out_col_corrected_data = tb.getcol('CORRECTED_DATA')
    tb.close()

    # Get input arrays
    tb.open(ms_in, nomodify=True)
    in_rows = tb.nrows()
    in_col_data = tb.getcol('DATA')
    in_col_model_data = tb.getcol('MODEL_DATA')
    in_col_corrected_data = tb.getcol('CORRECTED_DATA')
    in_col_ant1 = tb.getcol('ANTENNA1')
    in_col_ant2 = tb.getcol('ANTENNA2')
    in_col_weight = tb.getcol('WEIGHT')
    tb.close()

    num_antennas = _get_num_antennas(ms_out)
    num_baselines = num_antennas * (num_antennas - 1) / 2
    out_time_idx = numpy.zeros((num_baselines,), dtype=numpy.int32)

    for row in range(in_rows):
        a1 = in_col_ant1[row]
        a2 = in_col_ant2[row]
        weight = int(round(in_col_weight[0, row]))
        data = in_col_data[0, 0, row]
        model_data = in_col_model_data[0, 0, row]
        corrected_data = in_col_corrected_data[0, 0, row]
        b = a1 * (num_antennas - 1) - (a1 - 1) * a1 / 2 + a2 - a1 - 1;
        # print '%6i a1:%3i a2:%3i b:%6i w:%i:' % (row, a1, a2, b, weight),
        for t in range(out_time_idx[b], out_time_idx[b] + weight):
            out_row = t * num_baselines + b
            # print '(%2i)%6i' % (t, out_row),
            out_col_data[0, 0, out_row] = data
            out_col_model_data[0, 0, out_row] = model_data
            out_col_corrected_data[0, 0, out_row] = corrected_data
        # print ''
        out_time_idx[b] += weight

    tb.open(ms_out, nomodify=False)
    tb.putcol('DATA', out_col_data)
    tb.putcol('MODEL_DATA', out_col_model_data)
    tb.putcol('CORRECTED_DATA', out_col_corrected_data)
    tb.close()


if __name__ == '__main__':
    """
    iterate the data. and explode into copy of the reference for the given
    baseline for the number of time samples given by the weight in the bda ms.
    A updated start time index will have to be tracked for each baseline.
    """
    run()
