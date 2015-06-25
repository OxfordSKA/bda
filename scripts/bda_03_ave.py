#!/usr/bin/python
"""BDA with CASA mstransform task."""

import numpy as np
import os
import shutil
import time


def inv_sinc(arg):
    """Newton-Raphson method for calculating arcsinc(x), from Obit."""
    import numpy as np
    x1 = 0.001
    for i in range(0, 1000):
        x0 = x1
        a = x0 * np.pi
        x1 = x0 - ((np.sin(a) / a) - arg) / \
            ((a * np.cos(a) - np.pi * np.sin(a)) / (a**2))
        if (np.fabs(x1 - x0) < 1.0e-6):
            break
    return x1


def run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max):
    """Run the CASA mstransform task."""
    tb.open(ms_in + '/SPECTRAL_WINDOW', nomodify=True)
    freqs = tb.getcell('CHAN_FREQ', 0)
    tb.close()
    if len(freqs) != 1:
        print 'ERROR: can only use single channel data.'
        return
    freq = freqs[0]
    wavelength = 299792458.0 / freq
    delta_uv = inv_sinc(1.0 / max_fact) / (fov_radius * (np.pi / 180.))
    delta_uv *= wavelength  # convert to metres
    print ''
    print 'freq %.3e' % freq
    print 'fov  %.3f' % fov_radius
    print 'maxf %.3f' % max_fact
    print 'd_uv %.3e' % delta_uv
    print 'd_t  %s' % dt_max
    print ''
    if os.path.isdir(ms_out):
        shutil.rmtree(ms_out)
    mstransform(vis=ms_in,
                outputvis=ms_out,
                createmms=False,
                datacolumn='all',
                usewtspectrum=False,
                combinespws=False,
                chanaverage=False,
                hanning=False,
                timeaverage=True,
                timebin=dt_max,
                timespan='scan',
                maxuvwdistance=delta_uv)


if __name__ == "__main__":
    # -------------------------------------------------------------------------
    dt = 1.6  # Correlator dump time. TODO(BM) get this from the MS.
    dt_max = '%.2fs' % (5.0 * dt)  # Maximum allowed averaging time.
    max_fact = 1.01   # Maximum amplitude loss factor.
    fov_radius = 0.9  # Field of view radius (input into mstransform)
    # -------------------------------------------------------------------------

    # Average the model data.
    t0 = time.time()
    ms_in = os.path.join('vis', 'model.ms')
    ms_out = os.path.join('vis', 'model_mstransform.ms')
    if os.path.isdir(ms_in):
        run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max)
        print '+ Time taken in averaging = %.3fs [%s]' % \
            (time.time() - t0, ms_out)

    # Average the corrupted data.
    t0 = time.time()
    ms_in = os.path.join('vis', 'corrupted.ms')
    ms_out = os.path.join('vis', 'corrupted_ave.ms')
    if os.path.isdir(ms_in):
        run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max)
        print '+ Time taken in averaging = %.3fs [%s]' % \
            (time.time() - t0, ms_out)
