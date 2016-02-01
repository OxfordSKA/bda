# -*- coding: utf-8 -*-

import numpy
import os
from os.path import join
import shutil
import json
from bda import utilities


def _main():
    """Make a copy of the corrupted visibilities and add noise."""
    settings = utilities.byteify(json.load(open(config_file)))
    sim_dir = settings['path']
    ms_in = join(sim_dir, '%s.ms' % settings['ms_name']['corrupted'])
    if not os.path.isdir(ms_in):
        print 'ERROR: missing input visibility file. %s' % ms_in
        return
    ms_out = join(sim_dir, '%s_%s.ms' % (settings['ms_name']['corrupted'],
                                         settings['ms_modifier']['noisy']))
    if os.path.isdir(ms_out):
        return

    shutil.copytree(ms_in, ms_out)

    tb.open(ms_out, nomodify=True)
    num_rows = tb.nrows()
    col_data = tb.getcol('DATA')
    tb.close()

    numpy.random.seed(settings['noise']['seed'])
    noise_std = settings['noise']['noise_rms_jy']
    print 'Adding noise (%f Jy STD) to: %s...' % (noise_std, ms_out)
    noise = numpy.random.randn(num_rows) + 1j * numpy.random.randn(num_rows)
    noise *= noise_std

    col_data[0, 0, :] += noise

    tb.open(ms_out, nomodify=False)
    tb.putcol('DATA', col_data)
    tb.putcol('CORRECTED_DATA', col_data)
    tb.close()

    ms_out_2 = join(sim_dir, '%s_%s.ms' % (settings['ms_name']['model'],
                                           settings['ms_modifier']['reference']))
    print 'Adding noise (%f Jy STD) to: %s...' % (noise_std, ms_out_2)
    tb.open(ms_out_2, nomodify=False)
    tb.putcol('DATA', noise.reshape(1, 1, num_rows))
    tb.close()


if __name__ == '__main__':
    _main()
