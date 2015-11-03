# -*- coding: utf-8 -*-
"""Plot a summary of results from the BDA pipeline."""

import numpy
from os.path import join
import matplotlib.pyplot as plt
import pickle


def _load_gains(settings):
    sim_dir = settings['path']
    gains_file = '%s_%s.gains.pickle' % (settings['ms_name']['corrupted'],
                                         settings['ms_modifier']['sub_sampled'])
    gains_file = join(sim_dir, gains_file)
    gains = pickle.load(open(gains_file))
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
        ax.plot(x, numpy.abs(gains[i]),
                alpha=0.7, linewidth=0.5)
    ax.set_ylabel('Amplitude', fontsize='small')
    ax.grid(True)
    ax.set_ylim(1.0-0.03, 1.0+0.03)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax = axes2d[1]
    for i in range(len(gains)):
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


def _plot_corrupting_gains(settings):
    gains = _load_gains(settings)
    _plot_all_gains(gains, settings)
    antennas = _plot_gains(gains, settings)
    _plot_zoom_1(gains, settings, antennas[0],
                 'gains_zoom_ant_%i.png' % antennas[0])
    _plot_zoom_1(gains, settings, antennas[1],
                 'gains_zoom_ant_%i.png' % antennas[1])
    _plot_zoom_1(gains, settings, antennas[2],
                 'gains_zoom_ant_%i.png' % antennas[2])


def run(settings):
    _plot_corrupting_gains(settings)
