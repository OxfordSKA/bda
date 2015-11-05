# -*- coding: utf-8 -*-
"""Plot a summary of results from the BDA pipeline."""

import numpy
from os.path import join
import matplotlib.pyplot as plt
import pickle
from bda.utilities import adev
import copy


def _load_gains(settings):
    sim_dir = settings['path']
    gains_file = '%s_%s.gains.pickle' % (settings['ms_name']['corrupted'],
                                         settings['ms_modifier']['sub_sampled'])
    gains_file = join(sim_dir, gains_file)
    gains = pickle.load(open(gains_file))
    gains = copy.copy(gains)
    return gains


def _plot_all_gains(gains, settings):
    sim_dir = settings['path']
    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(6.5, 5.0))
    fig.subplots_adjust(left=0.10, bottom=0.10, right=0.97, top=0.95,
                        hspace=0.10, wspace=0.0)
    dt = (settings['sim']['observation']['dt_s'] /
          settings['sim']['observation']['over_sample'])
    x = numpy.arange(len(gains[0])) * dt

    ax = axes2d[0]
    for i in range(len(gains)):
        ax.plot(x, numpy.abs(gains[i]), alpha=0.7, linewidth=0.5)
    ax.set_ylabel('Amplitude', fontsize='small')
    ax.grid(True)
    ax.set_ylim(1.0-0.03, 1.0+0.03)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')

    ax = axes2d[1]
    for i in range(len(gains)):
        # print '%03i - %.3f' % (i,  numpy.std(gains[i]))
        ax.plot(x, numpy.degrees(numpy.angle(gains[i])),
                alpha=0.7, linewidth=0.5)
    ax.set_ylabel('Phase [degrees]', fontsize='small')
    ax.grid(True)
    ax.set_ylim(-160, 160)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.set_xlabel('Time [seconds]', fontsize='small')
    plt.savefig(join(sim_dir, 'gains_all.png'), dpi=300)


def _plot_gains(gains, settings):
    sim_dir = settings['path']
    # antennas = range(6)
    numpy.random.seed(20)
    antennas = numpy.random.randint(low=0, high=197, size=6)
    antennas = numpy.sort(antennas)
    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(6.5, 5.0))
    fig.subplots_adjust(left=0.1, bottom=0.10, right=0.97, top=0.95,
                        hspace=0.22, wspace=0.0)
    dt = (settings['sim']['observation']['dt_s'] /
          settings['sim']['observation']['over_sample'])
    x = numpy.arange(len(gains[0])) * dt
    ax = axes2d[0]
    for i in antennas:
        ax.plot(x, numpy.abs(gains[i]),
                linewidth=1.0,
                label='antenna %i' % i)
    ax.set_ylabel('Amplitude', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')

    ax = axes2d[1]
    for i in antennas:
        ax.plot(x, numpy.degrees(numpy.angle(gains[i])),
                linewidth=1.0,
                label='antenna %i' % i)
    ax.set_ylabel('Phase [degrees]', fontsize='small')
    ax.grid(True)
    # ax.set_ylim(-180, 180)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.set_xlabel('Time [seconds]', fontsize='small')

    ax = axes2d[1]
    ax.legend(bbox_to_anchor=(0, 1.02, 1.00, 0.5),
              loc=4,
              labelspacing=0.1,
              mode='expand',
              borderaxespad=0,
              handlelength=3.0,
              ncol=3,
              fontsize='x-small')

    plt.savefig(join(sim_dir, 'gains.png'), dpi=300)
    return antennas


def _plot_zoom_1(gains, settings, antenna, plot_name):
    sim_dir = settings['path']
    # ----------
    i0 = 0
    # i1 = i0 + settings['sim']['observation']['over_sample'] * 4
    i1 = len(gains[0])
    # ----------
    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(6.5, 5.0))
    fig.subplots_adjust(left=0.15, bottom=0.10, right=0.97, top=0.95,
                        hspace=0.10, wspace=0.0)
    dt = (settings['sim']['observation']['dt_s'] /
          settings['sim']['observation']['over_sample'])
    x = numpy.arange(i0, i1) * dt

    ax = axes2d[0]
    y = numpy.abs(gains[antenna][i0:i1])
    # y -= numpy.mean(y)
    y -= y[0]
    ax.plot(x, y, linewidth=0.5, marker='None',
            markersize=2.0, label='antenna %i' % antenna)
    ax.set_ylabel('$\Delta$ amplitude', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.set_title('Corrupting gains for antenna %i' % antenna, fontsize='small')

    ax = axes2d[1]
    y = numpy.degrees(numpy.angle(gains[antenna][i0:i1]))
    # y -= numpy.mean(y)
    y -= y[0]
    ax.plot(x, y, linewidth=0.5, marker='None', markersize=2.0,
            label='antenna %i' % antenna)
    ax.set_ylabel('$\Delta$ phase [degrees]', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.set_xlabel('Time [seconds]', fontsize='small')
    plt.savefig(join(sim_dir, plot_name), dpi=300)


def _plot_adev(gains, settings):
    sim_dir = settings['path']
    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(6.5, 5.0))
    fig.subplots_adjust(left=0.15, bottom=0.10, right=0.97, top=0.95,
                        hspace=0.22, wspace=0.0)
    dt = (settings['sim']['observation']['dt_s'] /
          settings['sim']['observation']['over_sample'])
    num_times = len(gains[0])
    x = numpy.arange(num_times) * dt
    antennas = range(1, len(gains))
    tau_ = numpy.arange(5, num_times / 3, 5) * dt
    shape = (len(antennas), tau_.shape[0])
    amp_allan_dev = numpy.empty(shape=shape, dtype='f8')
    amp_allan_dev_error = numpy.empty(shape=shape, dtype='f8')
    phase_allan_dev = numpy.empty(shape=shape, dtype='f8')
    phase_allan_dev_error = numpy.empty(shape=shape, dtype='f8')
    for ia, a in enumerate(antennas):
        for i, t in enumerate(tau_):
            adev_, adev_err_, _ = adev(numpy.abs(gains[a]), dt, t)
            amp_allan_dev[ia, i] = adev_
            amp_allan_dev_error[ia, i] = adev_err_
            adev_, adev_err_, _ = adev(numpy.angle(gains[a]) *
                                       (180.0 / numpy.pi), dt, t)
            phase_allan_dev[ia, i] = adev_
            phase_allan_dev_error[ia, i] = adev_err_

    mean_amp_allan_dev = numpy.mean(amp_allan_dev, axis=0)
    mean_phase_allan_dev = numpy.mean(phase_allan_dev, axis=0)

    # Original errors, using reported adev_error.
    #mean_amp_allan_dev_err = numpy.mean(amp_allan_dev_error, axis=0)
    #mean_phase_allan_dev_err = numpy.mean(phase_allan_dev_error, axis=0)

    # Errors using std.dev. of the Allan deviations.
    mean_amp_allan_dev_err = numpy.std(amp_allan_dev, axis=0)
    mean_phase_allan_dev_err = numpy.std(phase_allan_dev, axis=0)

    ax = axes2d[0]
    # for i, a in enumerate(antennas):
    #     ax.errorbar(tau_, amp_allan_dev[i, :], amp_allan_dev_error[i, :],
    #                 fmt='.-')
    ax.errorbar(tau_, mean_amp_allan_dev, mean_amp_allan_dev_err, fmt='.-')
    ax.set_xlabel('Tau [seconds]')
    ax.set_ylabel('Amplitude Allan dev.')
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    print ax.get_ylim()

    ax = axes2d[1]
    # for i, a in enumerate(antennas):
    #     ax.errorbar(tau_, phase_allan_dev[i, :],
    #                 phase_allan_dev_error[i, :], fmt='.-')
    ax.errorbar(tau_, mean_phase_allan_dev, mean_phase_allan_dev_err,
                fmt='.-')
    ax.set_xlabel('Tau [seconds]')
    ax.set_ylabel('Phase Allan dev. [deg]')
    plt.savefig(join(sim_dir, 'adev.png'), dpi=300)


def _plot_corrupting_gains(settings):
    gains = _load_gains(settings)
    _plot_all_gains(gains, settings)
    _plot_adev(gains, settings)
    antennas = _plot_gains(gains, settings)
    _plot_zoom_1(gains, settings, antennas[0],
                 'gains_zoom_ant_%i.png' % antennas[0])
    _plot_zoom_1(gains, settings, antennas[1],
                 'gains_zoom_ant_%i.png' % antennas[1])
    _plot_zoom_1(gains, settings, antennas[2],
                 'gains_zoom_ant_%i.png' % antennas[2])


def run(settings):
    _plot_corrupting_gains(settings)
