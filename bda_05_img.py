#!/usr/bin/python -u
# To ignore numpy errors:
#     pylint: disable=E1101
# To ignore undefined CASA classes
#     pylint: disable=undefined-variable

def fov_to_cellsize(fov, imsize):
    import numpy as np
    rmax = np.sin(np.array(fov, np.double)/2.*(np.pi/180.))
    inc = rmax / (0.5 * np.array(imsize))
    cell = np.arcsin(inc)*((180.*3600.)/np.pi)
    return cell.tolist()

def casa_image(ms, rootname, imsize, fov, ra0, dec0, wplanes=None):
    """Make an image using CASA
    http://casa.nrao.edu/docs/CasaRef/imager-Module.html#x636-6490002.5
    """
    import shutil
    import os

    if not os.path.isdir(os.path.dirname(rootname)):
        os.mkdir(os.path.dirname(rootname))

    cell = fov_to_cellsize(fov, imsize) # arcsec

    print '-'*80
    print '+ Size     : %i   x %i   [pixels]' % (imsize[0], imsize[1])
    print '+ FoV      : %.2f x %.2f [deg]' % (fov[0], fov[1])
    print '+ Cellsize : %.4f x %.4f [arcsec]' % (cell[0], cell[1])
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
    #psf = rootname+'_psf.img'
    im.makeimage(image=dirty, type='observed', verbose=False)
    #im.makeimage(image=psf, type='psf', verbose=False)
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
    ms = os.path.join('vis', 'test_cor_ave.ms')
    all_ms = [ms for ms in os.listdir('vis') if ms[-2:] == 'ms']
    rootname = 'images/cor_ave'
    img_size = [512, 512]
    img_fov = [0.3, 0.3] # deg
    # TODO(BM) load ra, dec from sky model
    ra0 = -90.35458487600000
    dec0 = -7.67112399060000
    wplanes = None
    # ---------------------------------------

    t0 = time.time()
    for ms in all_ms:
        rootname = os.path.join('images', ms[:-3])
        ms = os.path.join('vis', ms)
        print '+ Imaging with CASA ... [ms=%s -> %s]' % (ms, rootname)
        casa_image(ms, rootname, img_size, img_fov, ra0, dec0, wplanes)
        print '  - Finished imaging in %.3fs' % (time.time()-t0)

if __name__ == "__main__":
    """If running the CASA imager, run with:
            casa --nologger --nogui --log2term -c bda_05_img.py
    """
    main()
