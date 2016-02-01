#!venv/bin/python

import pyfits
import aplpy
import matplotlib.pyplot as plt
from os.path import join
import os.path
import numpy


def plot_scaling(target, weighting):
    root_dir = 'bda_results'
    obs_length = [5, 10, 30, 60, 120]
    types = {
        'default': '_',
        'noisy_default': '_noisy_',
        'bda': '_bda_',
        'noisy_bda': '_noisy_bda_',
        'expanded_bda': '_bda_expanded_bda_',
        'noisy_expanded_bda': '_noisy_bda_expanded_bda_'
    }
    results = {
        'default': numpy.zeros(len(obs_length)),
        'noisy_default': numpy.zeros(len(obs_length)),
        'bda': numpy.zeros(len(obs_length)),
        'noisy_bda': numpy.zeros(len(obs_length)),
        'expanded_bda': numpy.zeros(len(obs_length)),
        'noisy_expanded_bda': numpy.zeros(len(obs_length)),
        'noise': numpy.zeros(len(obs_length))
    }
    for i, t in enumerate(obs_length):
        sim_dir = 'SIM_%04is' % t

        for type_key in types:
            type = types[type_key]
            model_file = 'calibrated%sMODEL_DATA_%s_%s.fits' % \
                         (type, target, weighting)
            model_file = join(root_dir, sim_dir, model_file)
            calibrated_file = 'calibrated%sCORRECTED_DATA_%s_%s.fits' % \
                              (type, target, weighting)
            calibrated_file = join(root_dir, sim_dir, calibrated_file)
            print model_file, os.path.isfile(model_file)
            print calibrated_file, os.path.isfile(calibrated_file)

            model = pyfits.getdata(model_file)
            calibrated = pyfits.getdata(calibrated_file)
            diff = calibrated - model
            results[type_key][i] = numpy.std(diff)
        noise_file = 'model_ref_DATA_%s_%s.fits' % (target, weighting)
        noise_file = join(root_dir, sim_dir, noise_file)
        noise = pyfits.getdata(noise_file)
        results['noise'][i] = numpy.std(noise)
        print ''

    for k in results:
        print k, results[k]



    ax = plt.gca()
    # ax.set_yscale('log', nonposy='clip')
    # ax.set_xscale('log', nonposy='clip')
    # ax.plot(obs_length, results['noise'], '.-')
    ax.plot(obs_length, results['expanded_bda'], 'b-')
    ax.plot(obs_length, results['noisy_expanded_bda'], 'r-')
    # ax.plot(obs_length, results['default'], 'gx-')
    # ax.plot(obs_length, results['bda'], 'g-')



if __name__ == '__main__':
    fig = plt.figure(figsize=(6.5, 5.0))
    ax = fig.add_subplot(111)
    plot_scaling('1', 'n')
    plot_scaling('0', 'n')
    plt.show()
