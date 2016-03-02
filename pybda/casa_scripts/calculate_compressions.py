# -*- coding: utf-8 -*-# -*- coding: utf-8 -*-

import numpy
import os
import shutil
import time
import math
import subprocess

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

def run_bda(ms_in, ms_out, max_fact, fov_radius, dt_max):
    subprocess.call(["../src/bda_2", ms_in, ms_out,
        '%.5f' % max_fact, '%.3f' % fov_radius, '%i' % dt_max])

if __name__ == "__main__":

    ms = 'model.ms'
    #fh = open('mstransform_compression_log.txt', 'w')
    fh = open('bda_compression_log.txt', 'w')
    fov_radius = 0.9 # degrees
    dt_max = '15s'
    dt_max_val = 15
    f = numpy.array([1.001, 1.002, 1.003, 1.004, 1.005, 1.010, \
                     1.015, 1.020, 1.030, 1.050, 1.070, 1.100])

    for i in range(len(f)):
        max_fact = f[i]
        if os.path.isdir(ms):
            ms_out = 't_%.3f_%.1fdeg_%s.ms' % (max_fact, fov_radius, dt_max)
            t0 = time.time()
            #run_mstransform(ms, ms_out, max_fact, fov_radius, dt_max)
            run_bda(ms, ms_out, max_fact, fov_radius, dt_max_val)
            tb.open(ms, nomodify=True)
            rows_in = tb.nrows()
            tb.close()
            tb.open(ms_out, nomodify=True)
            rows_out = tb.nrows()
            tb.close()
            print '%i:%i -> 1:%.3f' % (rows_out, rows_in,
                                       rows_in/float(rows_out))
            fh.write('%.3f,%.3f\n' % (max_fact, rows_in/float(rows_out)))
            print '+ Time taken in averaging = %.3fs [%s]' % \
                (time.time() - t0, ms_out)

    fh.close()
