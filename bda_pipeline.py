#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function, absolute_import
import json
import os
import argparse
from os.path import join
from shutil import copyfile
from pybda import simple_simulator, imager
import matplotlib.pyplot as pyplot
from pybda.calibrate import calibrate
import numpy
import time


def main(config_file):
    try:
        config = json.load(open(config_file))
    except ValueError as e:
        print('Error: Failed to parse JSON config file.')
        print(e.message)
        return

    # Create the simulation directory, if required.
    sim_dir = config['path']
    if not os.path.isdir(sim_dir):
        os.makedirs(sim_dir)
    sim_config = join(sim_dir, os.path.basename(config_file))
    if not os.path.exists(sim_config):
        copyfile(config_file, sim_config)

    # Simulate visibilities
    vis = simple_simulator.simulate_2(config)

    # TODO-BM: BDA / expand (vis: model, data)

    # Calibrate the visibilities with StefCal
    calibrate(vis)

    # TODO-BM: BDA

    # Image visibilities
    # print('xx', len(config['imaging']['images']))
    model = imager.make_image(config, vis, image_id=1)
    dirty = imager.make_image(config, vis, image_id=2)
    corrected = imager.make_image(config, vis, image_id=3)
    diff_model_dirty = model - dirty
    diff_model_corrected = model - corrected

    fig = pyplot.figure(figsize=(20, 12))
    ax = fig.add_subplot(231)
    im1 = ax.imshow(model, interpolation='nearest')
    ax.figure.colorbar(im1, ax=ax)
    ax.set_title('model')

    ax = fig.add_subplot(232)
    im2 = ax.imshow(dirty, interpolation='nearest')
    ax.figure.colorbar(im2, ax=ax)
    ax.set_title('dirty')

    ax = fig.add_subplot(233)
    im3 = ax.imshow(corrected, interpolation='nearest')
    ax.figure.colorbar(im3, ax=ax)
    ax.set_title('corrected')

    ax = fig.add_subplot(234)
    im4 = ax.imshow(diff_model_dirty, interpolation='nearest')
    ax.figure.colorbar(im4, ax=ax)
    ax.set_title('model - dirty')

    ax = fig.add_subplot(235)
    im5 = ax.imshow(diff_model_corrected, interpolation='nearest')
    ax.figure.colorbar(im5, ax=ax)
    ax.set_title('model - corrected')

    pyplot.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='BDA pipeline main.',
                                     epilog='')
    parser.add_argument('config', type=str, nargs='?', help='JSON config file')
    args = parser.parse_args()
    if args.config is None:
        parser.print_usage()
        print('%s: error, too few arguments.' % os.path.basename(__file__))
        exit(1)
    if not os.path.isfile(args.config):
        print("Error: Config file '%s' not found!" % args.config)

    main(args.config)
