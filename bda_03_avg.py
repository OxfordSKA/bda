#!/usr/bin/python

import numpy as np
import os
import shutil

def inv_sinc(arg):
    """
    Newton-Raphson method for calculating arcsinc(x), from Obit.
    """
    import numpy as np
    x1 = 0.001
    for i in range(0, 1000):
        x0 = x1
        a = x0 * np.pi
        x1 = x0-((np.sin(a)/a)-arg) / ((a*np.cos(a) - np.pi*np.sin(a)) / (a**2))
        if (np.fabs(x1 - x0) < 1.0e-6): break
    return x1

def run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max):

    tb.open(ms_in+'/SPECTRAL_WINDOW', nomodify=True)
    freqs = tb.getcell('CHAN_FREQ', 0)
    tb.close()

    if len(freqs) != 1:
        print 'ERROR: can only use single channel data.'
        return

    freq = freqs[0]
    wavelength = 299792458.0 / freq
    delta_uv = inv_sinc(1.0 / max_fact) / (fov_radius * (np.pi/180.0))
    delta_uv *= wavelength

    if os.path.isdir(ms_out):
        shutil.rmtree(ms_out)

    mstransform(vis=ms_in,
                outputvis=ms_out,
                createmms=False,
                datacolumn='all',  # TODO(BM) check this value
                usewtspectrum=False,
                combinespws=False,
                chanaverage=False,
                hanning=False,
                timeaverage=True,
                timebin=dt_max,
                timespan='scan',
                maxuvwdistance=delta_uv)


def main():
    # ----------------------------------------
    ms_in = os.path.join('vis', 'test_cor.ms')
    ms_out = os.path.join('vis', 'test_cor_ave.ms')
    max_fact = 1.01   # Maximum amplitude factor by which a source can drop.
    dt = 1.6  # Correlator dump time. TODO(BM) get this from the MS.
    dt_max = '%.14fs' % (10.0*dt)  # Maximum allowed averaging time, as a
                                   # CASA string.
    fov_radius = 0.9  # Field of view radius (input into mstransform)
    # ----------------------------------------

    run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max)


if __name__ == "__main__":
    main()
