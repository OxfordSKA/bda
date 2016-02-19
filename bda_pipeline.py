#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function, absolute_import
import json
import os
import argparse
from os.path import join
from shutil import copyfile
from pybda import simple_simulator, imager, bda, calibrate
import matplotlib.pyplot as pyplot
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
    #calibrate.run_calibrate(vis)

    # Run BDA
    ave_data = bda.run_bda(config, vis, 'data')

    # Image visibilities
    # print('xx', len(config['imaging']['images']))
    model_images = imager.run_imager(config, vis, 'model')
    dirty_images = imager.run_imager(config, vis, 'data')
    #corrected = imager.run_imager(config, vis, 'corrected')
    #bda_model_images = imager.run_imager(config, ave_model, 'data')
    bda_dirty_images = imager.run_imager(config, ave_data, 'data')
    #bda_corrected = imager.run_imager(config, ave_model, 'corrected')
    diff_model_dirty = []
    for model, dirty in zip(model_images, dirty_images):
        diff_model_dirty.append(model - dirty)
    diff_data_bda_non_bda = []
    for dirty, bda_dirty in zip(dirty_images, bda_dirty_images):
        diff_data_bda_non_bda.append(dirty - bda_dirty)
    #diff_model_corrected = model - corrected

    # for i in range(len(model_images)):
    for i in range(1):

        fig = pyplot.figure(figsize=(16, 8))
        ax = fig.add_subplot(231)
        im1 = ax.imshow(model_images[i], interpolation='nearest')
        ax.figure.colorbar(im1, ax=ax)
        ax.set_title('model')

        ax = fig.add_subplot(232)
        im2 = ax.imshow(dirty_images[i], interpolation='nearest')
        ax.figure.colorbar(im2, ax=ax)
        ax.set_title('dirty')

        ax = fig.add_subplot(233)
        im3 = ax.imshow(bda_dirty_images[i], interpolation='nearest')
        ax.figure.colorbar(im3, ax=ax)
        ax.set_title('BDA')

        ax = fig.add_subplot(234)
        im4 = ax.imshow(diff_model_dirty[i], interpolation='nearest')
        ax.figure.colorbar(im4, ax=ax)
        ax.set_title('model - dirty')

        ax = fig.add_subplot(235)
        im5 = ax.imshow(diff_data_bda_non_bda[i], interpolation='nearest')
        ax.figure.colorbar(im5, ax=ax)
        ax.set_title('dirty - BDA dirty')

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
