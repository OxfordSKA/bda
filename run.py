#!venv/bin/python

import json
import os
import argparse
from os.path import isdir, join
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

    settings = json.load(open(args.config))
    run_steps = settings['steps']

    # Spawn a CASA process to work with and put the config_file variable to
    # the CASA variable list.
    casa = drivecasa.Casapy(working_dir=os.path.curdir,
                            casa_logfile=False,
                            echo_to_stdout=True)
    casa.run_script(["config_file = '{}'".format(args.config)])

    # Top level simulation output directory.
    sim_dir = settings['path']

    # Create a copy of the settings file.
    if not os.path.isdir(sim_dir):
        os.makedirs(sim_dir)
    copyfile(args.config, join(sim_dir, args.config))

    # Simulation.
    if "sim" in run_steps and not \
            isdir(join(sim_dir, settings['sim']['output_ms'])):
        simulate.run(settings)

    # Corrupt simulation.
    if "corrupt" in run_steps and not \
            isdir(join(sim_dir, settings['corrupt']['output_ms'])):
        casa.run_script_from_file('bda/casa_scripts/corrupt.py')

    # Baseline dependent averaging.
    if "bda" in run_steps and not \
            isdir(join(sim_dir, settings['baseline_average']['output_ms'][0])):
        casa.run_script_from_file('bda/casa_scripts/baseline_average.py')

    # Calibration.
    if "calibrate" in run_steps and not \
            isdir(join(sim_dir, settings['calibration']['output_ms'][0])):
        casa.run_script_from_file('bda/casa_scripts/calibration.py')

    # Image.
    if "image" in run_steps:
        casa.run_script_from_file('bda/casa_scripts/image.py')

    # Plot results.
    if "plot" in run_steps:
        plot.run(settings)

