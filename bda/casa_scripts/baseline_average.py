# -*- coding: utf-8 -*-# -*- coding: utf-8 -*-
"""BDA with CASA mstransform task."""

import numpy
import os
from os.path import join
import shutil
import time
import sys
import math
import json
from bda import utilities


def get_time_info(ms):
    """."""
    tb.open(ms, nomodify=True)
    times = numpy.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [numpy.min(times), numpy.max(times)]
    num_times = len(times)
    length = time_range[1] - time_range[0]
    dt = length / (num_times - 1)
    return num_times, time_range, length, dt


def inv_sinc(arg):
    """Newton-Raphson method for calculating arcsinc(x), from Obit."""
    x1 = 0.001
    for i in range(0, 1000):
        x0 = x1
        a = x0 * math.pi
        x1 = x0 - ((math.sin(a) / a) - arg) / \
            ((a * math.cos(a) - math.pi * math.sin(a)) / (a**2))
        if math.fabs(x1 - x0) < 1.0e-6:
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
    delta_uv = inv_sinc(1.0 / max_fact) / (fov_radius * (math.pi / 180.))
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

    settings = utilities.byteify(json.load(open(config_file)))
    # -------------------------------------------------------------------------
    sim_dir = settings['path']
    settings = settings['baseline_average']
    idt_max = settings['idt_max']
    max_fact = settings['max_fact']
    fov_radius = settings['fov_radius_deg']
    # -------------------------------------------------------------------------

    print '-' * 60
    print '+ Simulation directory:', sim_dir
    print '+ idt_max             :', idt_max
    print '+ max_fact            :', max_fact
    print '+ fov_radius          :', fov_radius
    print '-' * 60

    t_all = time.time()
    for ms in zip(settings['input_ms'], settings['output_ms']):
        t0 = time.time()
        ms_in = join(sim_dir, ms[0])
        ms_out = join(sim_dir, ms[1])
        if os.path.isdir(ms_in):
            _, _, _, dt = get_time_info(ms_in)
            dt_max = '%.5fs' % (idt_max * dt)
            print '-' * 60
            print '+ dt                  :', dt
            print '+ dt_max              :', dt_max
            print '-' * 60
            run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max)
            print '+ Time taken in averaging = %.3fs [%s]' % \
                (time.time() - t0, ms_out)

    print ''
    print '*' * 60
    print '+Total time taken = %.3fs' % (time.time() - t_all)
    print '*' * 60
