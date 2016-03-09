#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function, absolute_import
import json
import os
import argparse
from os.path import join
from shutil import copyfile
from pybda import simple_simulator, imager, bda, calibrate, expand_bda
#import matplotlib.pyplot as pyplot
import numpy
import time
import pyfits
import glob

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
    t0 = time.time()
    vis = simple_simulator.simulate_2(config)

    # Image model and corrupted data.
    original_model     = imager.run(config, vis, 'model', 'original_model')
    original_corrupted = imager.run(config, vis, 'data', 'original_corrupted')

    # Run initial BDA (compress) and overwrite the original data (expand).
    ave_model1, _ = bda.run(config, vis, 'model')
    expand_bda.run(vis['num_antennas'], ave_model1, 'data', vis, 'model')
    del ave_model1
    ave_data1, _ = bda.run(config, vis, 'data')
    expand_bda.run(vis['num_antennas'], ave_data1, 'data', vis, 'data')
    del ave_data1

    # Image corrupted data post-initial-BDA.
    post_bda1_corrupted = imager.run(config, vis, 'data', 'bda1_corrupted')

    # Calibrate the visibilities with StefCal.
    calibrate.run(vis, verbose=False)
    del vis['data']

    # Image calibrated data.
    calibrated = imager.run(config, vis, 'corrected', 'bda1_calibrated')

    # Run final BDA (compress).
    ave_model, _ = bda.run(config, vis, 'model')
    ave_data, compression_ratio = bda.run(config, vis, 'corrected')
    del vis

    # Image model and calibrated data post-final-BDA.
    post_bda2_model     = imager.run(config, ave_model, 'data', 'bda2_model')
    post_bda2_corrected = imager.run(config, ave_data, 'data', 'bda2_calibrated')
    del ave_model
    del ave_data

    # Save compression ratio.
    with open(join(sim_dir, 'compression.txt'), 'w') as f:
        f.write('%.4f\n' % (compression_ratio))

    # -------------------------------------------------------------------------
    # Make image differences
    num_images = len(config['imaging']['images'])
    diff_final_bda_corrected_final_bda_model = []
    diff_final_bda_model_original_model = []
    diff_original_corrupted_bda_corrupted = []
    for i in range(num_images):
        # final BDA corrected - final BDA model (for final metrics)
        diff = post_bda2_corrected[i] - post_bda2_model[i]
        diff_final_bda_corrected_final_bda_model.append(diff)
        fits_write(sim_dir, 'diff_corrected_and_model_'+str(i)+'.fits', diff)

        # final BDA model - original model (Bill Cotton's test)
        diff = post_bda2_model[i] - original_model[i]
        diff_final_bda_model_original_model.append(diff)
        fits_write(sim_dir, 'diff_models_'+str(i)+'.fits', diff)

        # original corrupted - first BDA/expand corrupted (see effect of BDA)
        diff = original_corrupted[i] - post_bda1_corrupted[i]
        diff_original_corrupted_bda_corrupted.append(diff)

    # Print time taken.
    print('- BDA pipeline took %.1f s' % (time.time() - t0))

    # for i in range(num_images):
    #     fig = pyplot.figure(figsize=(16, 11))
    #     ax = fig.add_subplot(3,3,1)
    #     im1 = ax.imshow(original_model[i], interpolation='nearest')
    #     ax.figure.colorbar(im1, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Original model')

    #     ax = fig.add_subplot(3,3,2)
    #     im2 = ax.imshow(original_corrupted[i], interpolation='nearest')
    #     ax.figure.colorbar(im2, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Original corrupted')

    #     ax = fig.add_subplot(3,3,3)
    #     im3 = ax.imshow(post_bda1_corrupted[i], interpolation='nearest')
    #     ax.figure.colorbar(im3, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Initial BDA corrupted')

    #     ax = fig.add_subplot(3,3,4)
    #     im4 = ax.imshow(post_bda2_model[i], interpolation='nearest')
    #     ax.figure.colorbar(im4, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Final BDA model')

    #     ax = fig.add_subplot(3,3,5)
    #     im5 = ax.imshow(post_bda2_corrected[i], interpolation='nearest')
    #     ax.figure.colorbar(im5, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Final BDA calibrated')

    #     ax = fig.add_subplot(3,3,6)
    #     im6 = ax.imshow(diff_final_bda_corrected_final_bda_model[i], interpolation='nearest')
    #     ax.figure.colorbar(im6, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Final BDA (calibrated - model)')

    #     ax = fig.add_subplot(3,3,7)
    #     im7 = ax.imshow(diff_final_bda_model_original_model[i], interpolation='nearest')
    #     ax.figure.colorbar(im7, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Model (Final BDA - Original) (Cotton)')

    #     ax = fig.add_subplot(3,3,8)
    #     im8 = ax.imshow(diff_original_corrupted_bda_corrupted[i], interpolation='nearest')
    #     ax.figure.colorbar(im8, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Corrupted (Original - Initial BDA)')

    #     ax = fig.add_subplot(3,3,9)
    #     im9 = ax.imshow(calibrated[i], interpolation='nearest')
    #     ax.figure.colorbar(im9, ax=ax)
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     ax.set_title('Calibrated')

    #     pyplot.show()


def fits_write(sim_dir, filename, image):
    pathname = join(sim_dir, filename)
    header = pyfits.header.Header([('SIMPLE', True), ('NAXIS', 2),
        ('NAXIS1', image.shape[0]), ('NAXIS2', image.shape[1])])
    if (os.path.exists(pathname)):
        os.remove(pathname)
    pyfits.writeto(pathname, image, header)


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
