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
import sys
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
    # print '(min = %.4f, max = %.4f mJy/beam)' % (np.min(d), np.max(d))
    if flux_scale == 1.e3:
        print '(min = %.4e, max = %.4e mJy/beam)' % (np.min(d), np.max(d))
    else:
        print '(min = %.4e, max = %.4e uJy/beam)' % (np.min(d), np.max(d))
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
    # f.axis_labels.hide_y()
    # f.tick_labels.hide_y()


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
    # f.axis_labels.hide_y()
    # f.tick_labels.hide_y()


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
    # f.axis_labels.hide_y()
    # f.tick_labels.hide_y()


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
    # f.axis_labels.hide_y()
    # f.tick_labels.hide_y()


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
    # f.axis_labels.hide_y()
    # f.tick_labels.hide_y()


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
    # f.axis_labels.hide_y()
    # f.tick_labels.hide_y()


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
    # f.axis_labels.hide_y()
    # f.tick_labels.hide_y()


def plot_single(hdu, title, stretch, c_label, cmin, cmax, fig,
                subplot=(1, 1, 1)):
    f = aplpy.FITSFigure(hdu, figure=fig, subplot=subplot)
    f.show_colorscale(vmin=cmin, vmax=cmax, cmap='seismic', stretch=stretch)
    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(c_label)
    f.colorbar.set_axis_label_font(size='small')
    f.colorbar.set_font(size='small')
    f.tick_labels.set_font(size='small')
    f.axis_labels.set_font(size='small')
    f.set_title(title, fontsize='medium')
    # 1 2 3
    # 4 5 6
    # 7 8
    if subplot[0] == 3 and subplot[1] == 3:
        if subplot[2] == 1:
            f.add_grid()
            f.grid.set_color('white')
            f.grid.set_alpha(0.3)
            f.axis_labels.hide_x()
            f.tick_labels.hide_x()
        if subplot[2] == 4:
            f.add_grid()
            f.grid.set_color('white')
            f.grid.set_alpha(0.3)
            f.axis_labels.hide_x()
            f.tick_labels.hide_x()
        if subplot[2] == 7:
            f.add_grid()
            f.grid.set_color('black')
            f.grid.set_alpha(0.3)
        if subplot[2] == 2:
            f.add_grid()
            f.grid.set_color('white')
            f.grid.set_alpha(0.3)
            f.axis_labels.hide_x()
            f.tick_labels.hide_x()
            f.axis_labels.hide_y()
            f.tick_labels.hide_y()
        if subplot[2] == 5:
            f.add_grid()
            f.grid.set_color('white')
            f.grid.set_alpha(0.3)
            f.axis_labels.hide_x()
            f.tick_labels.hide_x()
            f.axis_labels.hide_y()
            f.tick_labels.hide_y()
        if subplot[2] == 8:
            f.add_grid()
            f.grid.set_color('black')
            f.grid.set_alpha(0.3)
            f.axis_labels.hide_y()
            f.tick_labels.hide_y()
        if subplot[2] == 3:
            f.add_grid()
            f.grid.set_color('black')
            f.grid.set_alpha(0.3)
            f.axis_labels.hide_y()
            f.tick_labels.hide_y()
            f.axis_labels.hide_x()
            f.tick_labels.hide_x()
        if subplot[2] == 6:
            f.add_grid()
            f.grid.set_color('black')
            f.grid.set_alpha(0.3)
            f.axis_labels.hide_y()
            f.tick_labels.hide_y()
    if subplot[0] == 1 and subplot[1] == 2:
        if subplot[2] == 1:
            f.add_grid()
            f.grid.set_color('white')
            f.grid.set_alpha(0.3)
        if subplot[2] == 2:
            f.add_grid()
            f.grid.set_color('white')
            f.grid.set_alpha(0.3)
            f.axis_labels.hide_y()
            f.tick_labels.hide_y()



def main(sim_dir):
    """."""

    plot_gains = False
    plot_image_corrupted = False
    plot_image_grid = True

    if plot_image_grid:
        model = os.path.join(sim_dir, 'images', 'model_dirty.fits')
        model_bda = os.path.join(sim_dir, 'images', 'model_bda_dirty.fits')
        calibrated = os.path.join(sim_dir, 'images', 'calibrated_dirty.fits')
        calibrated_bda = os.path.join(sim_dir, 'images',
                                      'calibrated_bda_dirty.fits')

        # Plots:
        # 1 2 3
        # 4 5 6
        # 7 8
        fig = plt.figure(figsize=(20, 18))

        stretch = 'sqrt'
        c_label = r'mJy/beam'

        # 1st column : model
        cmin = -50.0  # mJy/beam
        cmax = 1000.0  # mJy/beam
        hdu_model = primary_hdu(model)
        hdu_model_bda = primary_hdu(model_bda)
        hdu_model_diff = primary_hdu_diff(model, model_bda, 1.e6)
        plot_single(hdu_model, 'model', stretch, c_label, cmin, cmax,
                    fig, (3, 3, 1))
        plot_single(hdu_model_bda, 'BDA model', stretch, c_label, cmin, cmax,
                    fig, (3, 3, 4))
        cmin = -30.   # uJy/beam
        cmax = 30.    # uJy/beam
        plot_single(hdu_model_diff, 'model - BDA model', 'linear', r'uJy/beam',
                    cmin, cmax,
                    fig, (3, 3, 7))
        print 'model diff: '
        print 'STD  : %e %s' % (np.std(hdu_model_diff.data), r'uJy/beam')
        print 'mean : %e %s' % (np.mean(hdu_model_diff.data), r'uJy/beam')
        print 'min  : %e %s' % (np.min(hdu_model_diff.data), r'uJy/beam')
        print 'max  : %e %s' % (np.max(hdu_model_diff.data), r'uJy/beam')


        # 2nd column : calibrated
        hdu_cal = primary_hdu(calibrated)
        hdu_cal_bda = primary_hdu(calibrated_bda)
        hdu_cal_diff = primary_hdu_diff(calibrated, calibrated_bda)
        cmin = -50.0  # mJy/beam
        cmax = 1000.0  # mJy/beam
        plot_single(hdu_cal, 'calibrated', stretch, c_label, cmin, cmax, fig,
                    (3, 3, 2))
        plot_single(hdu_cal_bda, 'BDA calibrated', stretch, c_label,
                    cmin, cmax, fig,
                    (3, 3, 5))
        cmin = -30.0  # mJy/beam
        cmax = 30.0   # mJy/beam
        # plot_338(hdu_cal_diff, 'calibrated - BDA calibrated',
        #          'linear', c_label, cmin, cmax, fig)
        plot_single(hdu_cal_diff, 'calibrated - BDA calibrated', 'linear',
                    c_label,
                    cmin, cmax, fig, (3, 3, 8))
        print 'calibrated diff: '
        print 'STD  : %e %s' % (np.std(hdu_cal_diff.data), r'mJy/beam')
        print 'mean : %e %s' % (np.mean(hdu_cal_diff.data), r'mJy/beam')
        print 'min  : %e %s' % (np.min(hdu_cal_diff.data), r'mJy/beam')
        print 'max  : %e %s' % (np.max(hdu_cal_diff.data), r'mJy/beam')



        # 3rd column: diffs
        cmin = -1.e-2  # uJy/beam
        cmax = 1.e-2   # uJy/beam
        hdu_diff = primary_hdu_diff(model, calibrated, 1.e6)
        hdu_diff_bda = primary_hdu_diff(model_bda, calibrated_bda)
        c_label = r'uJy/beam'
        plot_single(hdu_diff, 'model - calibrated', 'linear', c_label,
                    cmin, cmax, fig, (3, 3, 3))
        print 'model - calibrated diff: '
        print 'STD  : %e %s' % (np.std(hdu_diff.data), r'uJy/beam')
        print 'mean : %e %s' % (np.mean(hdu_diff.data), r'uJy/beam')
        print 'min  : %e %s' % (np.min(hdu_diff.data), r'uJy/beam')
        print 'max  : %e %s' % (np.max(hdu_diff.data), r'uJy/beam')

        c_label = r'mJy/beam'
        cmin = -30.0  # mJy/beam
        cmax = 30.0   # mJy/beam
        plot_single(hdu_diff_bda, 'BDA model - BDA calibrated', 'linear',
                    c_label,
                    cmin, cmax, fig, (3, 3, 6))
        print 'BDA model - BDA calibrated diff: '
        print 'STD  : %e %s' % (np.std(hdu_diff_bda.data), r'mJy/beam')
        print 'mean : %e %s' % (np.mean(hdu_diff_bda.data), r'mJy/beam')
        print 'min  : %e %s' % (np.min(hdu_diff_bda.data), r'mJy/beam')
        print 'max  : %e %s' % (np.max(hdu_diff_bda.data), r'mJy/beam')


        plt.tight_layout()
        plt.savefig(os.path.join(sim_dir, 'test_image_grid.png'))

    if plot_gains:
        print 'Plotting gains...'
        gains_pickle = os.path.join(sim_dir, 'vis', 'gains.pickle')
        g = pickle.load(open(gains_pickle))

        x = np.arange(0, 600) * 0.1
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,
                                       sharex=True, figsize=(20, 5))
        stations = range(0, 254)
        for s in stations:
            phase = np.angle(g[s]) * 180.0 / np.pi
            gains = np.abs(g[s])
            ax1.plot(x, phase, '-', lw=1.0)
            ax2.plot(x, gains, '-', lw=1.0)
        ax2.set_ylim([0.98, 1.02])
        ax1.set_ylim([-180, 180])
        ax2.set_xlabel('time [seconds]', fontsize='medium')
        ax2.set_ylabel('amplitude', fontsize='medium')
        ax1.set_ylabel('phase [degrees]', fontsize='medium')
        plt.tight_layout()
        plt.savefig(os.path.join(sim_dir, 'test_gains.png'))

        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,
                                       sharex=True, figsize=(20, 5))
        # rand_stations = np.random.randint(low=0, high=254, size=5)
        #
        # stations = [0]
        # stations.extend(rand_stations)
        stations = [0, 51, 91, 125, 135, 226]
        print stations
        for s in stations:
            phase = np.angle(g[s]) * 180.0 / np.pi
            gains = np.abs(g[s])
            ax1.plot(x, phase, '-', lw=1.0)
            ax2.plot(x, gains, '-', lw=1.0)
        ax2.set_ylim([0.98, 1.02])
        ax1.set_ylim([-180, 180])
        ax2.set_xlabel('time [seconds]', fontsize='medium')
        ax2.set_ylabel('amplitude', fontsize='medium')
        ax1.set_ylabel('phase [degrees]', fontsize='medium')
        plt.tight_layout()
        plt.savefig(os.path.join(sim_dir, 'test_gains_2.png'))

    if plot_image_corrupted:
        print '*' * 80
        print 'Plotting corrupted...'
        corrupted = os.path.join(sim_dir, 'images', 'corrupted_dirty.fits')
        model = os.path.join(sim_dir, 'images', 'model_dirty.fits')

        plt.figure(figsize=(20, 10))

        hdu_corrupted = primary_hdu(corrupted)
        hdu_model = primary_hdu(model)

        fig = plt.figure(figsize=(20, 10))

        plot_single(hdu_model, 'model', 'sqrt', 'mJy/beam',
                    -60.0, 1000.0, fig, (1, 2, 1))
        plot_single(hdu_corrupted, 'corrupted', 'sqrt', 'mJy/beam',
                    -60.0, 1000.0, fig, (1, 2, 2))

        plt.savefig(os.path.join(sim_dir, 'test_model_corrupted.png'))


if __name__ == "__main__":

    if len(sys.argv) - 1 < 1:
        print 'Usage:'
        print ('  $ python bda/sdp_f2f_plot.py '
               '<simulation dir>')
        sys.exit(1)

    sim_dir = sys.argv[-1]
    if not os.path.isdir(sim_dir):
        print 'ERROR: simulation directory not found!'
        sys.exit(1)

    print '-' * 60
    print 'Simulation directory:', sim_dir
    print '-' * 60

    main(sim_dir)
