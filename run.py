#!venv/bin/python
# -*- coding: utf-8 -*-

import json
import os
import argparse
from os.path import join
from shutil import copyfile
import drivecasa
from bda import simulate, plot


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
                            echo_to_stdout=True)
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

    # Average sub-sampled model and corrupted model.
    casa.run_script_from_file('bda/casa_scripts/average_ms.py')

    # TODO-BM add thermal noise.

    # Baseline dependent averaging.
    # casa.run_script_from_file('bda/casa_scripts/baseline_average.py')

    # Calibration.
    casa.run_script_from_file('bda/casa_scripts/calibration.py')

    # Image.
    casa.run_script_from_file('bda/casa_scripts/image.py')

    # Plot results.
    # plot.run(settings)

