# -*- coding: utf-8 -*-# -*- coding: utf-8 -*-
"""BDA with CASA mstransform task."""

import numpy as np
import os
import shutil
import time
import sys


def get_time_info(ms):
    """."""
    tb.open(ms, nomodify=True)
    times = np.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [np.min(times), np.max(times)]
    num_times = len(times)
    length = time_range[1] - time_range[0]
    dt = length / (num_times - 1)
    return num_times, time_range, length, dt


def inv_sinc(arg):
    """Newton-Raphson method for calculating arcsinc(x), from Obit."""
    x1 = 0.001
    for i in range(0, 1000):
        x0 = x1
        a = x0 * np.pi
        x1 = x0 - ((np.sin(a) / a) - arg) / \
            ((a * np.cos(a) - np.pi * np.sin(a)) / (a**2))
        if np.fabs(x1 - x0) < 1.0e-6:
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
    print '*' * 80
    print 'freq    %.3e' % freq
    print 'fov     %.3f' % fov_radius
    print 'maxf    %.3f' % max_fact
    print 'd_uv    %.3e' % delta_uv
    print 'd_t     %s' % dt_max
    print '*' * 80
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
    idt_max = 100
    max_fact = 1.01   # Maximum amplitude loss factor.
    fov_radius = 0.9  # Field of view radius (input into mstransform)
    # -------------------------------------------------------------------------

    if len(sys.argv) - 1 < 1:
        print 'Usage:'
        print ('  $ casa --nologger --nogui -c scripts/bda_02_cor.py '
               '<simulation dir>')
        sys.exit(1)

    sim_dir = sys.argv[-1]
    if not os.path.isdir(sim_dir):
        print 'ERROR: simulation directory not found!'
        sys.exit(1)

    print '-' * 60
    print '+ Simulation directory:', sim_dir
    print '+ idt_max             :', idt_max
    print '+ max_fact            :', max_fact
    print '+ fov_radius          :', fov_radius
    print '-' * 60

    # Average the model data.
    t_all = time.time()
    t0 = time.time()
    ms_in = os.path.join(sim_dir, 'vis', 'model.ms')
    ms_out = os.path.join(sim_dir, 'vis', 'model_mstransform.ms')
    _, _, _, dt = get_time_info(ms_in)
    dt_max = '%.2fs' % (idt_max * dt)
    print '-' * 60
    print '+ dt                  :', dt
    print '+ dt_max              :', dt_max
    print '-' * 60
    if os.path.isdir(ms_in):
        run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max)
        print '+ Time taken in averaging = %.3fs [%s]' % \
            (time.time() - t0, ms_out)

    # Average the corrupted data.
    t0 = time.time()
    ms_in = os.path.join(sim_dir, 'vis', 'corrupted.ms')
    ms_out = os.path.join(sim_dir, 'vis', 'corrupted_mstransform.ms')
    _, _, _, dt = get_time_info(ms_in)
    dt_max = '%.2fs' % (idt_max * dt)
    print '-' * 60
    print '+ dt                  :', dt
    print '+ dt_max              :', dt_max
    print '-' * 60
    if os.path.isdir(ms_in):
        run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max)
        print '+ Time taken in averaging = %.3fs [%s]' % \
            (time.time() - t0, ms_out)

    print ''
    print '*' * 60
    print '+Total time taken = %.3fs' % (time.time() - t_all)
    print '*' * 60
