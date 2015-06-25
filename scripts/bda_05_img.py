#!/usr/bin/python -u


def fov_to_cellsize(fov, imsize):
    import numpy as np
    rmax = np.sin(np.array(fov, np.double)/2.*(np.pi/180.))
    inc = rmax / (0.5 * np.array(imsize))
    cell = np.arcsin(inc)*((180.*3600.)/np.pi)
    return cell.tolist()


def casa_image(ms, rootname, data_column, imsize, fov, ra0, dec0, wplanes=None):
    """Make an image using CASA
    http://casa.nrao.edu/docs/CasaRef/imager-Module.html#x636-6490002.5
    """
    import shutil
    import os

    if not os.path.isdir(os.path.dirname(rootname)):
        os.mkdir(os.path.dirname(rootname))

    cell = fov_to_cellsize(fov, imsize) # arcsec

    print '-'*80
    print '+ Size     : %i pixels' % (imsize[0])
    print '+ FoV      : %.2f deg' % (fov[0])
    print '+ Cellsize : %.4f arcsec' % (cell[0])
    print '+ RA0      : %.4f deg' % (ra0)
    print '+ Dec0     : %.4f deg' % (dec0)
    print '-'*80

    im.open(ms, usescratch=False, compress=False)
    im.defineimage(nx=imsize[0], ny=imsize[1], cellx='%.12farcsec' % cell[0],
                   celly='%.12farcsec' % cell[1],
                   stokes='I', mode='mfs', step=1, spw=[-1], outframe='',
                   veltype='radio',
                   phasecenter=me.direction('J2000', '%.14fdeg' % ra0,
                                            '%.14fdeg' % dec0))
    # im.weight(type='natural')
    im.weight(type='uniform')
    if wplanes:
        im.setoptions(ftmachine='wproject', wprojplanes=wplanes,
                      gridfunction='SF', padding=1.2,
                      dopbgriddingcorrections=True, applypointingoffsets=False)
    else:
        im.setoptions(ftmachine='ft', gridfunction='SF', padding=1.2,
                      dopbgriddingcorrections=True, applypointingoffsets=False)

    dirty = rootname+'_dirty.img'
    # psf = rootname+'_psf.img'
    if data_column is 'DATA':
        im.makeimage(image=dirty, type='observed', verbose=False)  # DATA column
    elif data_column is 'CORRECTED_DATA':
        im.makeimage(image=dirty, type='corrected', verbose=False)  # CORRECTED_DATA column
    elif data_column is 'MODEL_DATA':
        im.makeimage(image=dirty, type='model', verbose=False)  # MODEL_DATA column
    else:
        print 'ERROR: Unknown data column!'
        return
    im.close()

    ia.open(dirty)
    ia.tofits(rootname+'_dirty.fits', overwrite=True)
    ia.close()
    # ia.open(psf)
    # ia.tofits(rootname+'_psf.fits', overwrite=True)
    # ia.close()
    if os.path.isdir(dirty):
        shutil.rmtree(dirty)
    # if os.path.isdir(psf):
    #     shutil.rmtree(psf)

def main():
    import os
    import time

    # ---------------------------------------
    # ms = os.path.join('vis', 'test_cor_ave.ms')
    all_ms = [ms for ms in os.listdir('vis') if ms[-2:] == 'ms']
    # rootname = 'images/cor_ave'
    img_size = [512, 512]
    img_fov = [0.01, 0.01]  # deg
    # TODO(BM) load ra, dec from sky model
    ra0 = -90.35458487600000
    dec0 = -7.67112399060000
    wplanes = None
    # ---------------------------------------


    for ms in all_ms:
        if not 'ave' in ms:
            continue
        rootname = os.path.join('images', ms[:-3])
        ms = os.path.join('vis', ms)

        if 'calibrated' in ms:
            data_column = 'CORRECTED_DATA'
        else:
            data_column = 'DATA'
        print '+ Imaging with CASA ... [ms=%s -> %s : %s]' % (ms, rootname,
                                                              data_column)
        t0 = time.time()
        casa_image(ms, rootname, data_column, img_size, img_fov, ra0, dec0,
                   wplanes)
        print '*'*80
        print '  - Finished imaging in %.3fs' % (time.time()-t0)
        print '*'*80


def main_temp():
    import os
    import time

    # ---------------------------------------
    img_size = [512, 512]
    img_fov = [0.01, 0.01]  # deg
    # TODO(BM) load ra, dec from sky model
    ra0 = -90.35458487600000
    dec0 = -7.67112399060000
    wplanes = None
    # ---------------------------------------

    # casa_image('vis/MS02.ms', 'images/MS02_DATA', 'DATA', 
    #         img_size, img_fov, ra0, dec0, wplanes)
    # casa_image('vis/MS02.ms', 'images/MS02_MODEL_DATA', 'MODEL_DATA', 
    #         img_size, img_fov, ra0, dec0, wplanes)
    # casa_image('vis/MS02.ms', 'images/MS02_CORRECTED_DATA', 'CORRECTED_DATA', 
    #         img_size, img_fov, ra0, dec0, wplanes)
    casa_image('vis/MS08.ms', 'images/MS08_DATA', 'DATA', 
            img_size, img_fov, ra0, dec0, wplanes)
    casa_image('vis/MS08.ms', 'images/MS08_MODEL_DATA', 'MODEL_DATA', 
            img_size, img_fov, ra0, dec0, wplanes)
    casa_image('vis/MS08.ms', 'images/MS08_CORRECTED_DATA', 'CORRECTED_DATA', 
            img_size, img_fov, ra0, dec0, wplanes)

if __name__ == "__main__":
    """If running the CASA imager, run with:
            casa --nologger --nogui --log2term -c bda_05_img.py
    """
    main()
    #main_temp()
