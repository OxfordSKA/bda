#!venv/bin/python
# -*- coding: utf-8 -*-

import json
import os
import argparse
from os.path import join
from shutil import copyfile
import drivecasa
from bda import simulate, plot_gains, plot_images
from bda.util.fits_diff import fits_diff


def _diff_cal_model(suffix, label, sim_dir):
    diff = 'diff%scal_model_%s.fits' % (label, suffix)
    cal = 'calibrated%sCORRECTED_DATA_%s.fits' % (label, suffix)
    model = 'calibrated%sMODEL_DATA_%s.fits' % (label, suffix)
    if not os.path.isfile(join(sim_dir, cal)):
        return
    fits_diff(join(sim_dir, diff), join(sim_dir, cal),
              join(sim_dir, model))


def _run(settings_file):
    try:
        settings = json.load(open(settings_file))
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
    casa.run_script(["config_file = '{}'".format(settings_file)])

    # Top level simulation output directory.
    sim_dir = settings['path']

    # Create the simulation directory.
    if not os.path.isdir(sim_dir):
        os.makedirs(sim_dir)

    # Create a copy of the settings file.
    if not os.path.exists(join(sim_dir, args.config)):
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
    # casa.run_script_from_file('bda/casa_scripts/add_noise.py')

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

    # Make some difference fits images.
    # TODO-BM move this into its own module.
    # weight 'u' = uniform, 'n' = natural
    # label: 0 = @ source, 1 = @ phase centre
    for label in ['0', '1']:
        for weight in ['u', 'n']:
            _suffix = '%s_%s' % (label, weight)
            for _label in ['_', '_noisy_', '_bda_', '_noisy_bda_',
                           '_bda_expanded_', '_noisy_bda_expanded_',
                           '_bda_expanded_bda_', '_noisy_bda_expanded_bda_']:
                _diff_cal_model(_suffix, _label, sim_dir)

    # Plot results fits images / diffs.
    plot_images.run(settings)


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

    _run(args.config)

