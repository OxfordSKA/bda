# -*- coding: utf-8 -*-

from pybda import utilities
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

    ms_names = ['%s_%s' % (settings['ms_name']['corrupted'],
                           settings['ms_modifier']['bda']),
                '%s_%s_%s' % (settings['ms_name']['corrupted'],
                              settings['ms_modifier']['noisy'],
                              settings['ms_modifier']['bda'])
                ]
    suffix = settings['ms_modifier']['expanded']

    ms_model = join(sim_dir, '%s.ms' % settings['ms_name']['model'])
    ms_ref = join(sim_dir, '%s_%s.ms' % (settings['ms_name']['model'],
                                         settings['ms_modifier']['reference']))

    tb.open(ms_model, nomodify=True)
    in_col_model_data_2 = tb.getcol('DATA')
    tb.close()

    for ms in ms_names:
        ms_in = join(sim_dir, '%s.ms' % ms)
        ms_out = join(sim_dir, '%s_%s.ms' % (ms, suffix))

        if os.path.isdir(ms_out) or not os.path.isdir(ms_in):
            continue

        print 'Expanding ms:', ms_in
        print '          ->:', ms_out

        # Copy reference
        shutil.copytree(ms_ref, ms_out)

        # Add model and corrected data columns to the measurement set.
        clearcal(vis=ms_out, addmodel=True)
        # Get output array(s) to write into.
        tb.open(ms_out, nomodify=True)
        out_col_uvw = tb.getcol('UVW')
        out_col_data = tb.getcol('DATA')
        out_col_model_data = tb.getcol('MODEL_DATA')
        out_col_corrected_data = tb.getcol('CORRECTED_DATA')
        tb.close()

        # Get input arrays
        tb.open(ms_in, nomodify=True)
        in_rows = tb.nrows()
        in_col_uvw = tb.getcol('UVW')
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
            weight = int(round(in_col_weight[:, row]))
            data = in_col_data[:, :, row]
            uvw = in_col_uvw[:, row]
            model_data = in_col_model_data[:, :, row]
            corrected_data = in_col_corrected_data[:, :, row]
            b = a1 * (num_antennas - 1) - (a1 - 1) * a1 / 2 + a2 - a1 - 1
            # print '%6i a1:%3i a2:%3i b:%6i w:%i:' % (row, a1, a2, b, weight),
            for t in range(out_time_idx[b], out_time_idx[b] + weight):
                out_row = t * num_baselines + b
                # print '(%2i)%6i' % (t, out_row),
                out_col_uvw[:, out_row] = uvw
                out_col_data[:, :, out_row] = data
                out_col_model_data[:, :, out_row] = model_data
                out_col_corrected_data[:, :, out_row] = corrected_data
            # print ''
            out_time_idx[b] += weight

        tb.open(ms_out, nomodify=False)
        # FIXME-BM: also copying uvw data values may be a bad idea...?
        # tb.putcol('UVW', out_col_uvw)
        tb.putcol('DATA', out_col_data)
        # NOTE: while it might be ok to calibrate against the smooth
        #       (non-bda-expanded model) we cannot difference against it
        #       without leaving stripes in the difference image.
        #       This is because the smooth model represents a different 'PSF'
        #       or 'forward simulated model' to the bda-expanded model.
        tb.putcol('MODEL_DATA', out_col_model_data)
        # tb.putcol('MODEL_DATA', in_col_model_data_2)
        tb.putcol('CORRECTED_DATA', out_col_corrected_data)
        tb.close()


if __name__ == '__main__':
    """
    Iterate the data and explode into copy of the reference for the given
    baseline for the number of time samples given by the weight in the bda ms.
    A updated start time index will have to be tracked for each baseline.
    """
    run()
