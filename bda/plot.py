# -*- coding: utf-8 -*-
"""Plot a summary of results from the BDA pipeline."""

import numpy
import os
from os.path import join
import matplotlib.pyplot as plt
import pyfits
import aplpy
import pickle


def adjust_header(header):
    for i, k in enumerate(header):
        if not k or k == 'COMMENT':
            header.pop(k)
    return header


def get_hdu(file_name, flux_scale=1.0e3):
    data, header = pyfits.getdata(file_name, header=True)
    header = adjust_header(header)
    data *= flux_scale
    return pyfits.PrimaryHDU(data, header), data


def get_hdu_diff(file1, file2, flux_scale=1.0e3):
    d1, h1 = pyfits.getdata(file1, header=True)
    d2, h2 = pyfits.getdata(file2, header=True)
    h1 = adjust_header(h1)
    assert d1.shape == d2.shape
    d = d1 - d2
    d *= flux_scale
    return pyfits.PrimaryHDU(d, h1), d

def plot_model_bda(settings):
    sim_dir = settings['path']
    plotting = settings['plotting']
    stretch = 'linear'
    fig = plt.figure(figsize=(6.5, 10.0))
    subplot = (2, 1, 1)
    hdu, data = get_hdu(join(sim_dir, plotting['model']))
    cmin = data.min()
    cmax = data.max()
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=subplot)
    f.add_label(0.05, 0.95, '(a)',
                relative=True,
                multialignment='center',
                fontsize=14,
                bbox=dict(color='white', alpha=0.75, edgecolor='white',
                          pad=5.0, linewidth=0.0))
    f.show_colorscale(vmin=cmin, vmax=cmax, stretch=stretch, cmap='afmhot')
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(r'mJy / beam')
    f.colorbar.set_axis_label_font(size='small')
    f.tick_labels.set_font(size='small')
    f.axis_labels.set_font(size='small')
    f.set_yaxis_coord_type('scalar')
    f.set_xaxis_coord_type('scalar')
    f.tick_labels.set_yformat('%.3f')
    f.tick_labels.set_xformat('%.3f')
    f.add_grid()
    f.grid.set_color('white')
    f.grid.set_linestyle('--')
    f.grid.set_alpha(0.3)
    f.set_title('model', fontsize='small', weight='bold')

    subplot=(2, 1, 2)
    hdu, data = get_hdu_diff(join(sim_dir, plotting['model_bda']),
                             join(sim_dir, plotting['model']),
                             1.0e6)
    cmin = data.min()
    cmax = data.max()
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=subplot)
    f.add_label(0.05, 0.95, '(b)',
                relative=True,
                multialignment='center',
                fontsize=14,
                bbox=dict(color='white', alpha=0.75, edgecolor='none',
                          pad=5.0, linewidth=0.0))
    f.show_colorscale(vmin=cmin, vmax=cmax, stretch=stretch, cmap='afmhot')
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(r'$\mathrm{\mu}$Jy / beam')
    f.colorbar.set_axis_label_font(size='small')
    f.tick_labels.set_font(size='small')
    f.axis_labels.set_font(size='small')
    f.set_yaxis_coord_type('scalar')
    f.set_xaxis_coord_type('scalar')
    f.tick_labels.set_yformat('%.3f')
    f.tick_labels.set_xformat('%.3f')
    f.add_grid()
    f.grid.set_color('white')
    f.grid.set_linestyle('--')
    f.grid.set_alpha(0.3)
    f.set_title('bda model - model', fontsize='small', weight='bold')
    plt.savefig(join(sim_dir, 'test.png'), dpi=300)
    plt.show()

def plot_model(settings):
    sim_dir = settings['path']
    plotting = settings['plotting']
    stretch = 'linear'
    fig = plt.figure(figsize=(6.5, 5.0))
    hdu, data = get_hdu(join(sim_dir, plotting['model']))
    cmin = data.min()
    cmax = data.max()
    f = aplpy.FITSFigure(hdu, figure=fig)
    f.show_colorscale(vmin=cmin, vmax=cmax, stretch=stretch, cmap='afmhot')
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(r'mJy / beam')
    f.colorbar.set_axis_label_font(size='small')
    f.tick_labels.set_font(size='x-small')
    f.axis_labels.set_font(size='x-small')
    # f.set_yaxis_coord_type('scalar')
    # f.set_xaxis_coord_type('scalar')
    # f.tick_labels.set_yformat('%.3f')
    # f.tick_labels.set_xformat('%.3f')
    f.add_grid()
    f.grid.set_color('white')
    f.grid.set_linestyle('--')
    f.grid.set_alpha(0.3)
    f.set_title('model', fontsize='small', weight='bold')
    plt.savefig(join(sim_dir, 'model.png'), dpi=300)


def plot_corrupting_gains(settings):
    sim_dir = settings['path']
    gains = pickle.load(open(join(sim_dir, settings['plotting']['gains'])))
    # Gains is a dict where the key is the antenna and the value is the list
    # of complex gains as a function of time.
    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(6.5, 5.0))
    fig.subplots_adjust(left=0.15, bottom=0.10, right=0.99, top=0.94,
                        hspace=0.1, wspace=0.0)
    ax = axes2d[0]
    for i in range(197):
        ax.plot(numpy.abs(gains[i]))
    ax.set_ylabel('Amplitude', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')

    ax = axes2d[1]
    for i in range(197):
        ax.plot(numpy.degrees(numpy.angle(gains[i])))
    ax.set_ylabel('Phase [degrees]', fontsize='small')
    ax.grid(True)
    ax.set_ylim(-90, 90)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    plt.savefig(join(sim_dir, 'gains_all.png'), dpi=300)

    antennas = range(5)
    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(6.5, 5.0))
    fig.subplots_adjust(left=0.15, bottom=0.10, right=0.99, top=0.85,
                        hspace=0.1, wspace=0.0)
    ax = axes2d[0]
    for i in antennas:
        ax.plot(numpy.abs(gains[i]), linewidth=0.5, label='antenna %i' % i)
    ax.set_ylabel('Amplitude', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.legend(bbox_to_anchor=(0, 1.02, 1.00, 0.5),
              loc=4,
              labelspacing=0.1,
              # mode='expand',
              borderaxespad=0,
              handlelength=3.0,
              ncol=3,
              fontsize='x-small')

    ax = axes2d[1]
    for i in antennas:
        ax.plot(numpy.degrees(numpy.angle(gains[i])), linewidth=0.5)
    ax.set_ylabel('Phase [degrees]', fontsize='small')
    ax.grid(True)
    # ax.set_ylim(-180, 180)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.set_xlabel('Time index', fontsize='small')
    plt.savefig(join(sim_dir, 'gains.png'), dpi=300)

    i0 = 100
    i1 = i0 + settings['sim']['observation']['over_sample'] * 4
    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(6.5, 5.0))
    fig.subplots_adjust(left=0.15, bottom=0.10, right=0.99, top=0.85,
                        hspace=0.1, wspace=0.0)
    ax = axes2d[0]
    for i in antennas:
        ax.plot(numpy.arange(i0, i1), numpy.abs(gains[i][i0:i1]),
                linewidth=0.5, marker='x',
                markersize=2.0, label='antenna %i' % i)
    ax.set_ylabel('Amplitude', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.legend(bbox_to_anchor=(0, 1.02, 1.00, 0.5),
              loc=4,
              labelspacing=0.1,
              # mode='expand',
              borderaxespad=0,
              handlelength=3.0,
              ncol=3,
              fontsize='x-small')

    ax = axes2d[1]
    for i in antennas:
        ax.plot(numpy.arange(i0, i1),
                numpy.degrees(numpy.angle(gains[i][i0:i1])), linewidth=0.5,
                marker='x', markersize=2.0)
    ax.set_ylabel('Phase [degrees]', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.set_xlabel('Time index', fontsize='small')
    plt.savefig(join(sim_dir, 'gains_zoom.png'), dpi=300)

    antennas = [4]
    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(6.5, 5.0))
    fig.subplots_adjust(left=0.15, bottom=0.10, right=0.99, top=0.85,
                        hspace=0.1, wspace=0.0)
    ax = axes2d[0]
    for i in antennas:
        ax.plot(numpy.arange(i0, i1), numpy.abs(gains[i][i0:i1]),
                linewidth=0.5, marker='x',
                markersize=2.0, label='antenna %i' % i)
    ax.set_ylabel('Amplitude', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.legend(bbox_to_anchor=(0, 1.02, 1.00, 0.5),
              loc=4,
              labelspacing=0.1,
              # mode='expand',
              borderaxespad=0,
              handlelength=3.0,
              ncol=3,
              fontsize='x-small')

    ax = axes2d[1]
    for i in antennas:
        ax.plot(numpy.arange(i0, i1),
                numpy.degrees(numpy.angle(gains[i][i0:i1])), linewidth=0.5,
                marker='x', markersize=2.0)
    ax.set_ylabel('Phase [degrees]', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='x-small')
    ax.tick_params(axis='both', which='major', labelsize='x-small')
    ax.set_xlabel('Time index', fontsize='small')
    plt.savefig(join(sim_dir, 'gains_zoom_2.png'), dpi=300)


def run(settings):
    # plot_model_bda(settings)
    plot_corrupting_gains(settings)
