#!venv/bin/python
# -*- coding: utf-8 -*-

import json
import os
import argparse
from os.path import join
from shutil import copyfile
import drivecasa
from bda import simulate, plot_gains
from bda.util import fits_diff

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='BDA script runner.',
                                     epilog='')
    parser.add_argument('config', type=str, nargs='?', help='JSON config file.')
    args = parser.parse_args()
    if args.config is None:
        parser.print_usage()
        print "%s: error: too few arguments" % os.path.basename(__file__)
        exit(1)
    if not os.path.isfile(args.config):
        print "Error: Config file '%s' not found!" % args.config
        exit(1)

    try:
        settings = json.load(open(args.config))
    except ValueError as e:
        print 'ERROR: FAILED TO PARSE JSON CONFIG FILE!!'
        print e.message
        exit(1)

    # Spawn a CASA process to work with.
    casa = drivecasa.Casapy(working_dir=os.path.curdir,
                            casa_logfile=False,
                            echo_to_stdout=True,
                            timeout=None)
    # Add config_file variable to the CASA variable list.
    casa.run_script(["config_file = '{}'".format(args.config)])

    # Top level simulation output directory.
    sim_dir = settings['path']

    # Create a copy of the settings file.
    if not os.path.isdir(sim_dir):
        os.makedirs(sim_dir)
    copyfile(args.config, join(sim_dir, args.config))

    # Simulation.
    simulate.run(settings, overwrite=False)

    # Introduce gain corruptions to the sub-sampled simulation.
    casa.run_script_from_file('bda/casa_scripts/corrupt.py')

    # Plot results.
    plot_gains.run(settings)

    # Average sub-sampled model and corrupted model to default interval.
    casa.run_script_from_file('bda/casa_scripts/average_ms.py')

    # Add thermal noise.
    # TODO-BM: for noise to work well we have to image a longer data set ...
    # the calibration error due to noise should be gaussian so should go away?
    # Otherwise have to play with the calibration interval in gaincal...?
    casa.run_script_from_file('bda/casa_scripts/add_noise.py')

    # Baseline dependent averaging.
    # TODO-BM different BDA schemes
    casa.run_script_from_file('bda/casa_scripts/baseline_average.py')

    # Expand the BDA MS back to to non-BDA sampling/
    casa.run_script_from_file('bda/casa_scripts/expand_bda_ms.py')

    # Calibration.
    casa.run_script_from_file('bda/casa_scripts/calibration.py')

    # Re-compress the BDA MS
    casa.run_script_from_file('bda/casa_scripts/baseline_average.py')

    # Image.
    casa.run_script_from_file('bda/casa_scripts/image.py')

    # # Make some difference fits images.
    weight = 'u'
    idx = 0
    fits_diff.fits_diff(join(sim_dir, 'diff_cal_model.fits'),
                        join(sim_dir, 'calibrated_'
                                      'CORRECTED_DATA_%i_%s.fits' % (idx, weight)),
                        join(sim_dir, 'calibrated_'
                                      'MODEL_DATA_%i_%s.fits' % (idx, weight)))

    fits_diff.fits_diff(join(sim_dir, 'diff_noisy_cal_model.fits'),
                        join(sim_dir, 'calibrated_noisy_'
                                      'CORRECTED_DATA_%i_%s.fits' % (idx, weight)),
                        join(sim_dir, 'calibrated_noisy_'
                                      'MODEL_DATA_%i_%s.fits' % (idx, weight)))

    fits_diff.fits_diff(join(sim_dir, 'diff_bda_cal_model.fits'),
                        join(sim_dir, 'calibrated_bda_'
                                      'CORRECTED_DATA_%i_%s.fits' % (idx, weight)),
                        join(sim_dir, 'calibrated_bda_'
                                      'MODEL_DATA_%i_%s.fits' % (idx, weight)))

    fits_diff.fits_diff(join(sim_dir, 'diff_noisy_bda_cal_model.fits'),
                        join(sim_dir, 'calibrated_noisy_bda_'
                                      'CORRECTED_DATA_%i_%s.fits' % (idx, weight)),
                        join(sim_dir, 'calibrated_noisy_bda_'
                                      'MODEL_DATA_%i_%s.fits' % (idx, weight)))

    fits_diff.fits_diff(join(sim_dir, 'diff_bda_expanded_cal_model.fits'),
                        join(sim_dir, 'calibrated_bda_expanded_'
                                      'CORRECTED_DATA_%i_%s.fits' % (idx, weight)),
                        join(sim_dir, 'calibrated_bda_expanded_'
                                      'MODEL_DATA_%i_%s.fits' % (idx, weight)))

    fits_diff.fits_diff(join(sim_dir, 'diff_noisy_bda_expanded_cal_model.fits'),
                        join(sim_dir, 'calibrated_noisy_bda_expanded_'
                                      'CORRECTED_DATA_%i_%s.fits' % (idx, weight)),
                        join(sim_dir, 'calibrated_noisy_bda_expanded_'
                                      'MODEL_DATA_%i_%s.fits' % (idx, weight)))

    fits_diff.fits_diff(join(sim_dir, 'diff_bda_expanded_bda_cal_model.fits'),
                        join(sim_dir, 'calibrated_bda_expanded_bda_'
                                      'CORRECTED_DATA_%i_%s.fits' % (idx, weight)),
                        join(sim_dir, 'calibrated_bda_expanded_bda_'
                                      'MODEL_DATA_%i_%s.fits' % (idx, weight)))

    fits_diff.fits_diff(join(sim_dir, 'diff_corrupted_noise.fits'),
                        join(sim_dir, 'corrupted_DATA_%i_%s.fits' % (idx, weight)),
                        join(sim_dir, 'corrupted_noisy_DATA_%i_%s.fits' % (idx, weight)))

    # TODO-BM: Plot results fits images / diffs.
    # plot.run(settings)
