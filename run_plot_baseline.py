#!venv/bin/python
# -*- coding: utf-8 -*-

import os
import argparse
from os.path import join
from shutil import copyfile
import drivecasa


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot visibilities per '
                                                 'baseline.',
                                     epilog='')
    parser.add_argument('sim_dir',
                        type=str, nargs='?', help='Simulation directory.')
    args = parser.parse_args()
    if args.sim_dir is None:
        parser.print_usage()
        print "%s: error: too few arguments" % os.path.basename(__file__)
        exit(1)

    # Spawn a CASA process to work with.
    casa = drivecasa.Casapy(working_dir=os.path.curdir,
                            casa_logfile=False,
                            echo_to_stdout=True)

    # Introduce gain corruptions to the sub-sampled simulation.
    casa.run_script(["sim_dir = '{}'".format(args.sim_dir)])
    casa.run_script_from_file('bda/casa_scripts/plot_baseline.py')

