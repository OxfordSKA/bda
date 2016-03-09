#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function, absolute_import
import json
import os
import argparse
from os.path import join
import numpy
import pyfits
import glob


def get_stats(sim_dir):
    # Load the config file.
    config_file = glob.glob(join(sim_dir, '*.json'))
    try:
        config = json.load(open(config_file[0]))
    except ValueError as e:
        print('Error: Failed to parse JSON config file.')
        print(e.message)
        return

    # Get configuration values.
    tau                = config['corrupt']['tau_s']
    hurst_amp          = config['corrupt']['amplitude']['hurst']
    adev_amp           = config['corrupt']['amplitude']['allan_dev']
    hurst_phase        = config['corrupt']['phase']['hurst']
    adev_phase         = config['corrupt']['phase']['allan_dev']
    num_times          = config['sim']['observation']['num_times']
    max_fact           = config['baseline_average']['max_fact']
    fov_radius_deg     = config['baseline_average']['fov_radius_deg']
    max_average_time_s = config['baseline_average']['max_average_time_s']
    compression_ratio  = numpy.loadtxt(join(sim_dir, 'compression.txt'))

    # Get stats.
    diff_files = glob.glob(join(sim_dir, 'diff_*.fits'))
    for i in range(len(diff_files)):
        diff_name = diff_files[i]
        image = pyfits.getdata(diff_name)
        rms_diff = numpy.sqrt(numpy.mean(numpy.square(image)))
        min_diff = numpy.min(image)
        max_diff = numpy.max(image)
        with open(os.path.splitext(diff_name)[0] + '.txt', 'w') as f:
            f.write('%.4f, %.2f, %.1f, %.3f, '
                '%i, %.1f, %.1f, %.1e, %.1f, %.1e, '
                '%.6f, %.6f, %.6f\n' % (
                max_fact, fov_radius_deg, max_average_time_s, compression_ratio,
                num_times, tau, hurst_amp, adev_amp, hurst_phase, adev_phase,
                rms_diff, min_diff, max_diff))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get stats for BDA pipeline.',
                                     epilog='')
    parser.add_argument('sim_dir', type=str, nargs='?', help='Simulation dir')
    args = parser.parse_args()
    if args.sim_dir is None:
        parser.print_usage()
        print('%s: error, too few arguments.' % os.path.basename(__file__))
        exit(1)
    if not os.path.isdir(args.sim_dir):
        print("Error: Simulation dir '%s' not found!" % args.sim_dir)

    get_stats(args.sim_dir)
