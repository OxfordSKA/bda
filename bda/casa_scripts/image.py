# -*- coding: utf-8 -*-
"""Make images using CASA."""

import numpy
import math
import shutil
import os
import time
from os.path import join
import json
from bda import utilities

def fov_to_cellsize(fov, im_size):
    """Obatin cellsize from fov and image size."""
    r_max = numpy.sin(numpy.array(fov, numpy.double) / 2. * (numpy.pi / 180.))
    inc = r_max / (0.5 * numpy.array(im_size))
    cell = numpy.arcsin(inc) * ((180. * 3600.) / numpy.pi)
    return cell.tolist()


def casa_image(ms, rootname, data_column, imsize, fov, ra0, dec0,
               weighting, w_planes=None):
    """Make an image using CASA.

    http://casa.nrao.edu/docs/CasaRef/imager-Module.html#x636-6490002.5
    """
    if not os.path.isdir(os.path.dirname(rootname)):
        os.mkdir(os.path.dirname(rootname))

    cell = fov_to_cellsize(fov, imsize)  # arcsec

    print '-' * 80
    print '+ Size     : %i pixels' % (imsize[0])
    print '+ FoV      : %.2f deg' % (fov[0])
    print '+ Cellsize : %.4f arcsec' % (cell[0])
    print '+ RA0      : %.4f deg' % (ra0)
    print '+ Dec0     : %.4f deg' % (dec0)
    print '-' * 80

    im.open(ms, usescratch=False, compress=False)
    im.defineimage(nx=imsize[0], ny=imsize[1], cellx='%.12farcsec' % cell[0],
                   celly='%.12farcsec' % cell[1],
                   stokes='I', mode='mfs', step=1, spw=[-1], outframe='',
                   veltype='radio',
                   phasecenter=me.direction('J2000', '%.14fdeg' % ra0,
                                            '%.14fdeg' % dec0))
    # im.weight(type='natural')
    im.weight(type=weighting)
    if w_planes:
        im.setoptions(ftmachine='wproject', wprojplanes=w_planes,
                      gridfunction='SF', padding=1.2,
                      dopbgriddingcorrections=True, applypointingoffsets=False)
    else:
        im.setoptions(ftmachine='ft', gridfunction='SF', padding=1.2,
                      dopbgriddingcorrections=True, applypointingoffsets=False)

    dirty = rootname + '_dirty.img'
    # psf = rootname + '_psf.img'
    if data_column == 'DATA':
        # DATA column
        im.makeimage(image=dirty, type='observed', verbose=False)
    elif data_column == 'CORRECTED_DATA':
        # CORRECTED_DATA column
        im.makeimage(image=dirty, type='corrected', verbose=False)
    elif data_column == 'MODEL_DATA':
        # MODEL_DATA column
        im.makeimage(image=dirty, type='model', verbose=False)
    else:
        print 'ERROR: Unknown data column!'
        return
    im.close()
    ia.open(dirty)
    ia.tofits(rootname + '.fits', overwrite=True)
    ia.close()
    # ia.open(psf)
    # ia.tofits(rootname+'_psf.fits', overwrite=True)
    # ia.close()
    if os.path.isdir(dirty):
        shutil.rmtree(dirty)
    # if os.path.isdir(psf):
    #     shutil.rmtree(psf)


def _run():
    settings = utilities.byteify(json.load(open(config_file)))

    if 'imaging' not in settings:
        return

    sim_dir = settings['path']
    ms_files = [f for f in os.listdir(os.path.abspath(sim_dir))
                if f.endswith('.ms') and os.path.isdir(join(sim_dir, f))]

    settings = settings['imaging']
    column_spec = settings['columns']
    for f in ms_files:
        ms = join(sim_dir, f)
        for k in column_spec.keys():
            if k in ms:
                columns = column_spec[k]
                for column in columns:
                    root_name = os.path.splitext(ms)[0] + '_{}'.format(column)
                    t0 = time.time()
                    for i, image in enumerate(settings['images']):
                        image_name = '%s_%i' % (root_name, i)
                        if image['weighting'] == 'natural':
                            image_name += '_n'
                        elif image['weighting'] == 'uniform':
                            image_name += '_u'
                        if not os.path.exists(image_name + '.fits'):
                            print '+ Imaging: %s [%s] (%i)' % (ms, column, i)
                            casa_image(ms, image_name,
                                       column, image['size'], image['fov_deg'],
                                       image['ra_deg'], image['dec_deg'],
                                       image['weighting'], image['w_planes'])
                            print '*' * 80
                            print '  - Finished imaging in %.3fs' % \
                                  (time.time() - t0)
                            print '*' * 80

if __name__ == '__main__':
    _run()
