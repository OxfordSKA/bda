# -*- coding: utf-8 -*-
"""Plot a summary of results from the BDA pipeline."""


import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyfits
import aplpy
import pickle
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.colors import LogNorm


def adjust_header(header):
    """Adjust a FITS header for use with aplpy."""
    for i, k in enumerate(header):
        if not k or k == 'COMMENT':
            header.pop(k)
    return header


def primary_hdu(file_name, flux_scale=1.e3):
    """."""
    data, header = pyfits.getdata(file_name, header=True)
    header = adjust_header(header)
    data *= flux_scale
    print 'Loading %s ' % file_name
    print '(min = %.4f, max = %.4f mJy/beam)' % (np.min(data), np.max(data))
    return pyfits.PrimaryHDU(data, header)


def primary_hdu_diff(file_1, file_2, flux_scale=1.e3):
    """."""
    d1, h1 = pyfits.getdata(file_1, header=True)
    d2, h2 = pyfits.getdata(file_2, header=True)
    h1 = adjust_header(h1)
    assert(d1.shape == d2.shape)
    d = d1 - d2
    d *= flux_scale
    print 'Loading %s - %s' % (file_1, file_2)
    print '(min = %.4f, max = %.4f mJy/beam)' % (np.min(d), np.max(d))
    print '(min = %.4e, max = %.4e mJy/beam)' % (np.min(d), np.max(d))
    return pyfits.PrimaryHDU(d, h1)


def plot_331(hdu, title, stretch, c_label, cmin, cmax, fig):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=(3, 3, 1))
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')


def plot_332(hdu, title, stretch, c_label, cmin, cmax, fig):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=(3, 3, 2))
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')
    f.axis_labels.hide_y()
    f.tick_labels.hide_y()


def plot_333(hdu, title, stretch, c_label, cmin, cmax, fig):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=(3, 3, 3))
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')
    f.axis_labels.hide_y()
    f.tick_labels.hide_y()


def plot_334(hdu, title, stretch, c_label, cmin, cmax, fig):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=(3, 3, 4))
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')
    f.axis_labels.hide_y()
    f.tick_labels.hide_y()


def plot_335(hdu, title, stretch, c_label, cmin, cmax, fig):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=(3, 3, 5))
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')
    f.axis_labels.hide_y()
    f.tick_labels.hide_y()


def plot_336(hdu, title, stretch, c_label, cmin, cmax, fig):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=(3, 3, 6))
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')
    f.axis_labels.hide_y()
    f.tick_labels.hide_y()


def plot_337(hdu, title, stretch, c_label, cmin, cmax, fig):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=(3, 3, 7))
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')
    f.axis_labels.hide_y()
    f.tick_labels.hide_y()


def plot_338(hdu, title, stretch, c_label, cmin, cmax, fig):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=(3, 3, 8))
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')
    f.axis_labels.hide_y()
    f.tick_labels.hide_y()

def plot_single(hdu, title, stretch, c_label, cmin, cmax, fig,
                subplot=(1, 1, 1)):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=subplot)
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='xx-small')
    f.colorbar.set_font(size='xx-small')
    f.tick_labels.set_font(size='xx-small')
    f.axis_labels.set_font(size='small')
    # f.add_grid()
    f.set_title(title, fontsize='small')



def main(sim_dir):
    """."""

    plot_image_grid = False
    plot_image_corrupted = True
    plot_gains = False

    # ave_type = 'bda'
    # sim_dir = 'out_new'
    ave_type = 'mstransform'
    sim_dir = 'out_default'

    if plot_image_grid:
        model = os.path.join(sim_dir, 'images', 'model_dirty.fits')
        model_bda = os.path.join(sim_dir, 'images',
                                 'model_%s_dirty.fits' % ave_type)
        calibrated = os.path.join(sim_dir, 'images', 'calibrated_dirty.fits')
        calibrated_bda = os.path.join(sim_dir, 'images',
                                      'calibrated_%s_dirty.fits' % ave_type)

        # Plots:
        # 1 2 3
        # 4 5 6
        # 7 8
        fig = plt.figure(figsize=(20, 20))

        stretch = 'linear'
        c_label = r'mJy/beam'

        # 1st column : model
        cmin = -50.0  # mJy/beam
        cmax = 500.0  # mJy/beam
        hdu_model = primary_hdu(model)
        hdu_model_bda = primary_hdu(model_bda)
        hdu_model_diff = primary_hdu_diff(model, model_bda)
        plot_331(hdu_model, 'model', stretch, c_label, cmin, cmax, fig)
        plot_334(hdu_model_bda, 'BDA model', stretch, c_label, cmin, cmax, fig)
        # cmin = -0.066  # mJy/beam
        # cmax = 0.11    # mJy/beam
        cmin = -0.066  # mJy/beam
        cmax = 0.11    # mJy/beam
        plot_337(hdu_model_diff, 'model - BDA model',
                 stretch, c_label, cmin, cmax, fig)

        # 2nd column : calibrated
        hdu_cal = primary_hdu(calibrated)
        hdu_cal_bda = primary_hdu(calibrated_bda)
        hdu_cal_diff = primary_hdu_diff(calibrated, calibrated_bda)
        cmin = -50.0  # mJy/beam
        cmax = 500.0  # mJy/beam
        plot_332(hdu_cal, 'calibrated', stretch, c_label, cmin, cmax, fig)
        plot_335(hdu_cal_bda, 'BDA calibrated', stretch, c_label, cmin, cmax, fig)
        cmin = -10.0  # mJy/beam
        cmax = 20.0   # mJy/beam
        plot_338(hdu_cal_diff, 'calibrated - BDA calibrated',
                 stretch, c_label, cmin, cmax, fig)

        # 3rd column: diffs
        cmin = -1.e-2  # uJy/beam
        cmax = 1.e-2   # uJy/beam
        hdu_diff = primary_hdu_diff(model, calibrated, 1.e6)
        hdu_diff_bda = primary_hdu_diff(model_bda, calibrated_bda)
        c_label = r'uJy/beam'
        plot_333(hdu_diff, 'model - calibrated',
                 stretch, c_label, cmin, cmax, fig)
        c_label = r'mJy/beam'
        cmin = -10.0  # mJy/beam
        cmax = 20.0   # mJy/beam
        plot_336(hdu_diff_bda, 'BDA model - BDA calibrated',
                 stretch, c_label, cmin, cmax, fig)

        plt.tight_layout()
        plt.savefig('test_image_grid.png')

    if plot_gains:
        print 'Plotting gains...'
        sim_dir = 'f2f'
        gains_pickle = os.path.join(sim_dir, 'vis', 'gains.pickle')
        g = pickle.load(open(gains_pickle))

        x = np.arange(0, 1000) * 0.1
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,
                                       sharex=True, figsize=(20, 5))
        stations = range(0, 254)
        # np.random.seed(5)
        # stations = np.random.randint(low=0, high=254, size=5)
        # print stations
        for s in stations:
            phase = np.angle(g[s]) * 180.0 / np.pi
            gains = np.abs(g[s])
            ax1.plot(x, phase, '-', lw=1.0)
            ax2.plot(x, gains, '-', lw=1.0)
        ax2.set_ylim([0.97, 1.03])
        ax1.set_ylim([-180, 180])
        ax2.set_xlabel('time [seconds]')
        ax2.set_ylabel('amplitude')
        ax1.set_ylabel('phase [degrees]')
        plt.tight_layout()
        plt.savefig('test_gains.png')

    if plot_image_corrupted:
        corrupted = os.path.join(sim_dir, 'images', 'corrupted_dirty.fits')
        model = os.path.join(sim_dir, 'images', 'model_dirty.fits')

        plt.figure(figsize=(20, 10))

        hdu_corrupted = primary_hdu(corrupted)
        hdu_model = primary_hdu(model)

        fig = plt.figure(figsize=(20, 10))

        plot_single(hdu_model, 'model', 'linear', 'mJy/beam',
                    -80.0, 300.0, fig, (1, 2, 1))
        plot_single(hdu_corrupted, 'corrupted', 'linear', 'mJy/beam',
                    -80.0, 300.0, fig, (1, 2, 2))

        plt.savefig('test_model_corrupted.png')


if __name__ == "__main__":
    sim_dir = 'out_default'
    main(sim_dir)
