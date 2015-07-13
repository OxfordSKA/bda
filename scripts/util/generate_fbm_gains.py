# -*- coding: utf-8 -*-
"""Generate fractional Brownian motion gain table for MATLAB."""


import sys
import scipy.io
import os
import numpy as np
from collections import OrderedDict
sys.path.append(os.path.join(os.getcwd(), 'scripts'))
from bda_00_util import eval_complex_gain


def main():
    """Generate gain table for MATLAB"""
    # ----------------------------------
    num_antennas = 254
    num_times = 800
    dt = 0.02  # seconds
    tau = round(1.0 / dt) * dt  # ~1.0 seconds.
    amp_std_t0 = 0.005  # Initial (t=0) gain amplitude variation per antenna
    amp_hurst = 0.8
    amp_allan_dev_fbm = 1.0e-4 / 4.0  # Allan deviation of gain amp @ tau
    amp_sigma_wn = 0.0
    phase_std_t0 = 45.0  # Initial (t=0) gain phase variation per antenna
    phase_hurst = 0.75
    phase_allan_dev_fbm = 2.0 / 4.0  # Allan deviation of gain phase @ tau
    phase_sigma_wn = 0.0
    # np.random.seed(666)
    gains_file = 'fbm_gains.mat'
    # ----------------------------------
    print 'Tau = %.2f' % tau

    gains = np.ones((num_times, num_antennas), dtype='c16')
    for i in range(1, num_antennas):
        g = eval_complex_gain(num_times, dt, amp_hurst, amp_allan_dev_fbm,
                              amp_sigma_wn, phase_hurst, phase_allan_dev_fbm,
                              phase_sigma_wn, amp_std_t0, phase_std_t0, tau)
        gains[:, i] = g

    gain_table = OrderedDict()
    gain_table['fbm_gains'] = {
        'gains': gains,
        'dt': dt,
        'num_times': num_times,
        'num_antennas': num_antennas,
        'tau': tau,
        'amp_std_t0': amp_std_t0,
        'amp_allen_dev_fbm': amp_allan_dev_fbm,
        'phase_std_t0': phase_std_t0,
        'phase_hurst': phase_hurst,
        'phase_allen_dev_fbm': phase_allan_dev_fbm
    }
    scipy.io.savemat(gains_file, gain_table)

    # Verify in MATLAB with:
    #
    # >> load('fbm_gains.mat')
    # >> plot(1:800, abs(fbm_gains.gains(:, 1:10)))
    # >> plot(1:800, angle(fbm_gains.gains(:, 1:10)))
    #

if __name__ == '__main__':
    main()


