# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import
import time
import numpy
from oskar._bda_utils import expand


def run(num_antennas, vis_compressed, input_name, vis_original, output_name):
    print('- Expanding compressed data...')
    t0 = time.time()
    expand(num_antennas, vis_compressed, input_name, vis_original, output_name)

    # num_baselines = num_antennas * (num_antennas - 1) / 2
    # num_input_vis = len(vis_compressed[input_name])
    # print('  - No. input visibilities  : %i' % num_input_vis)
    # out_time_idx = numpy.zeros((num_baselines,), dtype=numpy.int32)
    # for row in range(num_input_vis):
    #     a1 = vis_compressed['antenna1'][row]
    #     a2 = vis_compressed['antenna2'][row]
    #     weight = int(round(vis_compressed['weight'][row]))
    #     data = vis_compressed[input_name][row]
    #     b = a1 * (num_antennas - 1) - (a1 - 1) * a1 / 2 + a2 - a1 - 1
    #     for t in range(out_time_idx[b], out_time_idx[b] + weight):
    #         out_row = t * num_baselines + b
    #         vis_original[output_name][out_row] = data
    #     out_time_idx[b] += weight

    print('  - Visibilities expanded in %.2f s' % (time.time() - t0))
