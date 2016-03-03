#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function, absolute_import
import json
import os
import argparse
from os.path import join
from shutil import copyfile
from pybda import simple_simulator, imager, bda, calibrate, expand_bda
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

    # Simulate visibilities.
    vis = simple_simulator.simulate_2(config)

    # Image model and corrupted data.
    original_model     = imager.run_imager(config, vis, 'model')
    original_corrupted = imager.run_imager(config, vis, 'data')

    # Run initial BDA (compress) and overwrite the original data (expand).
    ave_model1 = bda.run_bda(config, vis, 'model')
    ave_data1 = bda.run_bda(config, vis, 'data')
    expand_bda.run_expand_bda(vis['num_antennas'], 
        ave_model1, 'data', vis, 'model')
    expand_bda.run_expand_bda(vis['num_antennas'], 
        ave_data1, 'data', vis, 'data')

    # Image corrupted data post-initial-BDA.
    post_bda1_corrupted = imager.run_imager(config, vis, 'data')

    # Calibrate the visibilities with StefCal.
    calibrate.run_calibrate(vis, verbose=False)

    # Image calibrated data.
    calibrated = imager.run_imager(config, vis, 'corrected')

    # Run final BDA (compress).
    ave_model = bda.run_bda(config, vis, 'model')
    ave_data = bda.run_bda(config, vis, 'corrected')

    # Image model and calibrated data post-final-BDA.
    post_bda2_model     = imager.run_imager(config, ave_model, 'data')
    post_bda2_corrected = imager.run_imager(config, ave_data, 'data')

    # Make image differences
    # -------------------------------------------------------------------------
    # final BDA corrected - final BDA model (for final metrics)
    diff_final_bda_corrected_final_bda_model = []
    for corr, model in zip(post_bda2_corrected, post_bda2_model):
        diff_final_bda_corrected_final_bda_model.append(corr - model)

    # final BDA model - original model (Bill Cotton's test)
    diff_final_bda_model_original_model = []
    for model_bda, model_original in zip(post_bda2_model, original_model):
        diff_final_bda_model_original_model.append(model_bda - model_original)

    # original corrupted - first BDA/expand corrupted (see effect of BDA)
    diff_original_corrupted_bda_corrupted = []
    for corr_orig, corr_bda in zip(original_corrupted, post_bda1_corrupted):
        diff_original_corrupted_bda_corrupted.append(corr_orig - corr_bda)

    # for i in range(len(original_model)):
    for i in range(1):
        fig = pyplot.figure(figsize=(16, 11))
        ax = fig.add_subplot(3,3,1)
        im1 = ax.imshow(original_model[i], interpolation='nearest')
        ax.figure.colorbar(im1, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Original model')

        ax = fig.add_subplot(3,3,2)
        im2 = ax.imshow(original_corrupted[i], interpolation='nearest')
        ax.figure.colorbar(im2, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Original corrupted')

        ax = fig.add_subplot(3,3,3)
        im3 = ax.imshow(post_bda1_corrupted[i], interpolation='nearest')
        ax.figure.colorbar(im3, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Initial BDA corrupted')

        ax = fig.add_subplot(3,3,4)
        im4 = ax.imshow(post_bda2_model[i], interpolation='nearest')
        ax.figure.colorbar(im4, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Final BDA model')

        ax = fig.add_subplot(3,3,5)
        im5 = ax.imshow(post_bda2_corrected[i], interpolation='nearest')
        ax.figure.colorbar(im5, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Final BDA calibrated')

        ax = fig.add_subplot(3,3,6)
        im6 = ax.imshow(diff_final_bda_corrected_final_bda_model[i], interpolation='nearest')
        ax.figure.colorbar(im6, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Final BDA (calibrated - model)')

        ax = fig.add_subplot(3,3,7)
        im7 = ax.imshow(diff_final_bda_model_original_model[i], interpolation='nearest')
        ax.figure.colorbar(im7, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Model (Final BDA - Original) (Cotton)')

        ax = fig.add_subplot(3,3,8)
        im8 = ax.imshow(diff_original_corrupted_bda_corrupted[i], interpolation='nearest')
        ax.figure.colorbar(im8, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Corrupted (Original - Initial BDA)')

        ax = fig.add_subplot(3,3,9)
        im9 = ax.imshow(calibrated[i], interpolation='nearest')
        ax.figure.colorbar(im9, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Calibrated')

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
