# -*- coding: utf-8 -*-# -*- coding: utf-8 -*-
"""BDA with CASA mstransform task."""

import numpy
import os
from os.path import join
import shutil
import time
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
    print 'd_uv    %.3e m' % delta_uv
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
    bda_settings = settings['baseline_average']
    idt_max = bda_settings['idt_max']
    max_fact = bda_settings['max_fact']
    fov_radius = bda_settings['fov_radius_deg']
    # -------------------------------------------------------------------------

    fh = open(join(sim_dir, 'bda_log.txt'))

    # model.ms -> model_bda.ms
    # corrupted.ms -> corrupted_bda.ms
    # corrupted_noisy.ms -> corrupted_noisy_bda.ms

    # calibrated_bda_expanded.ms -> calibrated_bda_expanded_bda.ms
    # calibrated_noisy_bda_expanded -> calibrated_noisy_bda_expanded_bda.ms

    ms_names = ['%s' % settings['ms_name']['model'],
                '%s' % settings['ms_name']['corrupted'],
                '%s_%s' % (settings['ms_name']['corrupted'],
                           settings['ms_modifier']['noisy']),
                '%s_%s_%s' % (settings['ms_name']['calibrated'],
                              settings['ms_modifier']['bda'],
                              settings['ms_modifier']['expanded']),
                '%s_%s_%s_%s' % (settings['ms_name']['calibrated'],
                                 settings['ms_modifier']['noisy'],
                                 settings['ms_modifier']['bda'],
                                 settings['ms_modifier']['expanded']),
                ]
    bda_suffix = settings['ms_modifier']['bda']

    for ms in ms_names:
        ms_in = join(sim_dir, '%s.ms' % ms)
        ms_out = join(sim_dir, '%s_%s.ms' % (ms, bda_suffix))
        if os.path.isdir(ms_out):
            continue
        if os.path.isdir(ms_in):
            _, _, _, dt = get_time_info(ms_in)
            dt_max = '%.5fs' % (idt_max * dt)
            print '-' * 60
            print '+ %s' % ms_in
            print '+ %s' % ms_out
            print '+ dt                  :', dt
            print '+ dt_max              :', dt_max
            print '-' * 60
            t0 = time.time()
            run_mstransform(ms_in, ms_out, max_fact, fov_radius, dt_max)
            tb.open(ms_in, nomodify=True)
            rows_in = tb.nrows()
            tb.close()
            tb.open(ms_out, nomodify=True)
            rows_out = tb.nrows()
            tb.close()
            print '%i:%i -> 1:%.3f' % (rows_out, rows_in,
                                       rows_in/float(rows_out))
            fh.write('%s,%s,%.3f' % (ms_in, ms_out, rows_in/float(rows_out)))
            print '+ Time taken in averaging = %.3fs [%s]' % \
                (time.time() - t0, ms_out)

    fh.close()
