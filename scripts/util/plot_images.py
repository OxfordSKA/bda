#!/usr/bin/python

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import pyfits
import aplpy

def apl_plot(filename):
    fig = aplpy.FITSFigure(filename)
    # fig.show_colorscale()
    fig.show_colorscale(cmap='seismic')
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text('Flux (Jy/beam')
    fig.colorbar.set_axis_label_font(size=12)
    # fig.add_grid()
    # fig.grid.set_xspacing(0.02)  # degrees
    # fig.grid.set_yspacing(0.02)  # degrees
    # fig.add_beam()  # needs BMAJ, BMIN?
    # fig.beam.show()
    # fig.set_theme('publication')
    fig.set_theme('pretty')
    fig.save(filename[:-5]+'.eps')
    # fig.save(filename[:-5]+'.png', transparent=True)

def db_image_plot(ax, data, title):
    data = 10.0*np.log10(np.abs(data)/np.nanmax(np.abs(data)))
    plt.imshow(data, interpolation='nearest', cmap='seismic')
    plt.clim([-30, 0])
    plt.title(title, fontsize=10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='4%', pad=0.05)
    cbar = plt.colorbar(cax=cax)
    plt.tick_params(axis='both', which='major', labelsize=10)
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])

def image_plot(ax, data, title):
    plt.imshow(data, interpolation='nearest', cmap='seismic')
    plt.title(title, fontsize=10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='4%', pad=0.05)
    cbar = plt.colorbar(cax=cax)
    plt.tick_params(axis='both', which='major', labelsize=10)
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])



def main():
    # ------------------------------------------
    model = os.path.join('images', 'model_dirty.fits')
    corrupted = os.path.join('images', 'corrupted_dirty.fits')
    calibrated = os.path.join('images', 'calibrated_dirty.fits')
    ave_model = os.path.join('images', 'model_ave_dirty.fits')
    ave_corrupted = os.path.join('images', 'corrupted_ave_dirty.fits')
    ave_calibrated = os.path.join('images', 'calibrated_ave_dirty.fits')
    # ------------------------------------------

    hdulist = pyfits.open(model)
    model_img = np.squeeze(hdulist[0].data)
    hdulist.close()
    hdulist = pyfits.open(corrupted)
    corrupted_img = np.squeeze(hdulist[0].data)
    hdulist.close()
    hdulist = pyfits.open(calibrated)
    calibrated_img = np.squeeze(hdulist[0].data)
    hdulist.close()
    hdulist = pyfits.open(ave_model)
    ave_model_img = np.squeeze(hdulist[0].data)
    hdulist.close()
    hdulist = pyfits.open(ave_corrupted)
    ave_corrupted_img = np.squeeze(hdulist[0].data)
    hdulist.close()
    hdulist = pyfits.open(ave_calibrated)
    ave_calibrated_img = np.squeeze(hdulist[0].data)
    hdulist.close()


    # apl_plot(model)
    # apl_plot(corrupted)

    fig = plt.figure(figsize=(20, 10))

    # Not averaged
    ax = fig.add_subplot(341, aspect='equal')
    db_image_plot(ax, model_img, 'Model')
    ax = fig.add_subplot(342, aspect='equal')
    db_image_plot(ax, corrupted_img, 'Corrupted')
    ax = fig.add_subplot(343, aspect='equal')
    db_image_plot(ax, calibrated_img, 'Calibrated')
    ax = fig.add_subplot(344, aspect='equal')
    diff = model_img-calibrated_img
    image_plot(ax, diff, 'Model-Calibrated (%.1e +/- %.1e)' %
               (np.mean(diff), np.std(diff)))
    cmin = np.mean(diff)-3.0*np.std(diff)
    cmax = np.mean(diff)+3.0*np.std(diff)
    plt.clim(cmin, cmax)
    # Averaged
    ax = fig.add_subplot(345, aspect='equal')
    db_image_plot(ax, ave_model_img, 'Ave: model')
    ax = fig.add_subplot(346, aspect='equal')
    db_image_plot(ax, ave_corrupted_img, 'Ave: corrupted')
    ax = fig.add_subplot(347, aspect='equal')
    db_image_plot(ax, ave_calibrated_img, 'Ave: calibrated')
    ax = fig.add_subplot(348, aspect='equal')
    diff = ave_model_img-ave_calibrated_img
    image_plot(ax, diff, 'Ave: Model-Calibrated (%.1e +/- %.1e)' %
               (np.mean(diff), np.std(diff)))
    cmin = np.mean(diff)-3.0*np.std(diff)
    cmax = np.mean(diff)+3.0*np.std(diff)
    plt.clim(cmin, cmax)

    ax = fig.add_subplot(349, aspect='equal')
    diff = model_img - ave_model_img
    image_plot(ax, diff, 'Model-ave_model (%.1e +/- %.1e)' %
               (np.mean(diff), np.std(diff)))
    cmin = np.mean(diff)-3.0*np.std(diff)
    cmax = np.mean(diff)+3.0*np.std(diff)
    plt.clim(cmin, cmax)

    ax = fig.add_subplot(3, 4, 10, aspect='equal')
    diff = corrupted_img - ave_corrupted_img
    image_plot(ax, diff, 'Corrupted-ave_corrupted (%.1e +/- %.1e)' %
               (np.mean(diff), np.std(diff)))
    cmin = np.mean(diff)-3.0*np.std(diff)
    cmax = np.mean(diff)+3.0*np.std(diff)
    plt.clim(cmin, cmax)

    ax = fig.add_subplot(3, 4, 11, aspect='equal')
    diff = calibrated_img - ave_calibrated_img
    image_plot(ax, diff, 'calibrated-ave_calibrated (%.1e +/- %.1e)' %
               (np.mean(diff), np.std(diff)))
    cmin = np.mean(diff)-3.0*np.std(diff)
    cmax = np.mean(diff)+3.0*np.std(diff)
    plt.clim(cmin, cmax)

    ax = fig.add_subplot(3, 4, 12, aspect='equal')
    ax.axis('off')
    plt.text(0.0, 1.0,
        "No averaging:\n"
        "  max(model) = %.3e\n"
        "  max(corrupted) = %.3e\n"
        "  max(calibrated) = %.3e\n"
        "  std(model-calibrated) = %.3e\n"
        "  max(model)-max(calibrated) = %.3e\n"
        "\n"
        "Baseline averaged:\n"
        "  max(model) = %.3e\n"
        "  max(corrupted) = %.3e\n"
        "  max(calibrated) = %.3e\n"
        "  std(model-calibrated) = %.3e\n"
        "  max(model)-max(calibrated) = %.3e\n"
        % (
           np.nanmax(model_img),
           np.nanmax(corrupted_img),
           np.nanmax(calibrated_img),
           np.nanstd(model_img-calibrated_img),
           (np.nanmax(model_img)-np.nanmax(calibrated_img)),
           np.nanmax(ave_model_img),
           np.nanmax(ave_corrupted_img),
           np.nanmax(ave_calibrated_img),
           np.nanstd(ave_model_img-ave_calibrated_img),
           (np.nanmax(ave_model_img)-np.nanmax(ave_calibrated_img)),
        ),
        size=10, ha='left', va='top')

    plt.savefig(os.path.join('images', 'bda_sim_images.eps'))
    plt.savefig(os.path.join('images', 'bda_sim_images.png'))
    plt.show()


if __name__ == "__main__":
    main()
