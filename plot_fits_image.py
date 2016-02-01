#!venv/bin/python
# -*- coding: utf-8 -*-

import aplpy
import matplotlib.pyplot as plt
import pyfits
import sys


def _adjust_header(header):
    for i, k in enumerate(header):
        if not k or k == 'COMMENT' or k == 'HISTORY':
            header.pop(k)
    return header


def _get_hdu(file_name, flux_scale=1.0e3):
    data, header = pyfits.getdata(file_name, header=True)
    header = _adjust_header(header)
    data *= flux_scale
    return pyfits.PrimaryHDU(data, header), data


if __name__ == '__main__':
    if len(sys.argv) - 1 != 1:
        print 'Usage: plot_fits_image.py <fits file>'
        sys.exit(1)

    fits_file = sys.argv[1]

    hdu, data = _get_hdu(fits_file, flux_scale=1.0e3)
    cmin = data.min()
    cmax = data.max()
    cmax = max(abs(cmin), abs(cmax))
    cmin = -cmax
    print cmin, cmax

    fig = plt.figure(figsize=(6.5, 6.5))
    stretch = 'linear'
    subplot = (1, 1, 1)

    f = aplpy.FITSFigure(hdu, figure=fig, subplot=subplot)
    # f.show_colorscale(vmin=cmin, vmax=cmax, stretch=stretch,
    #                   cmap='afmhot')
    f.show_colorscale(vmin=cmin, vmax=cmax, stretch=stretch,
                      cmap='seismic')
    # f.show_colorscale(vmax=1000.0, stretch=stretch, cmap='seismic')
    # f.show_colorscale(cmap='seismic')


    f.add_colorbar()
    f.colorbar.set_width(0.1)
    f.colorbar.set_axis_label_text(ur'mJy / beam')
    f.colorbar.set_axis_label_font(size='small')
    f.tick_labels.set_font(size='small')
    f.axis_labels.set_font(size='small')
    # f.set_yaxis_coord_type('scalar')
    # f.set_xaxis_coord_type('scalar')
    # f.tick_labels.set_yformat('%.3f')
    # f.tick_labels.set_xformat('%.3f')
    f.add_grid()
    f.grid.set_color('white')
    f.grid.set_linestyle('--')
    f.grid.set_alpha(0.3)
    f.set_title(fits_file, fontsize='small', weight='bold')
    f.save(fits_file + '.png', dpi=300)

    plt.show()
