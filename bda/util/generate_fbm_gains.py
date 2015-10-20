# -*- coding: utf-8 -*-
"""Generate fractional Brownian motion gain table for MATLAB."""


import scipy.io
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt


def main():
    """Generate gain table for MATLAB"""
    # ----------------------------------
    plot_gains = True

    num_antennas = 254
    num_times = 800
    dt = 0.1  # seconds
    tau = round(1.0 / dt) * dt  # ~1.0 seconds.

    amp_std_t0 = 0.003  # Initial (t=0) gain amplitude variation per antenna
    amp_hurst = 0.5
    amp_allan_dev_fbm = 1.e-3  # Allan deviation of gain amp @ tau
    amp_sigma_wn = 0.0

    phase_std_t0 = 30.0  # Initial (t=0) gain phase variation per antenna
    phase_hurst = 0.5
    phase_allan_dev_fbm = 0.5  # Allan deviation of gain phase @ tau
    phase_sigma_wn = 0.0

    smooth_gains = True
    window_length = 5
    window = 'hanning'

    np.random.seed(666)
    gains_file = 'fbm_gains_%i_%is.mat' % (num_times, dt * num_times)
    # ----------------------------------
    print 'Tau = %.2f' % tau

    gains = np.ones((num_times, num_antennas), dtype='c16')
    for i in range(1, num_antennas):
        g = eval_complex_gain(num_times, dt, amp_hurst, amp_allan_dev_fbm,
                              amp_sigma_wn, phase_hurst, phase_allan_dev_fbm,
                              phase_sigma_wn, amp_std_t0, phase_std_t0, tau)
        if smooth_gains:
            g = smooth(g, window_len=window_length, window=window)
        gains[:, i] = g[0:num_times]

    gain_table = OrderedDict()
    gain_table['fbm_gains_%i_%is' % (num_times, int(dt * num_times))] = {
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

    if plot_gains:
        antennas = range(1, 10)

        t = np.arange(0, num_times) * dt
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(221)
        ax.plot(t, np.abs(gains[:, antennas]))
        ax.set_ylabel('amplitude')
        ax.set_xlabel('Time [seconds]')
        ax = fig.add_subplot(222)
        ax.plot(t, np.angle(gains[:, antennas]) * (180. / np.pi))
        ax.set_ylim([-180, 180])
        ax.set_ylabel('phase [degrees]')
        ax.set_xlabel('Time [seconds]')

        tau_ = np.arange(1, num_times / 3, 5) * dt
        shape = (len(antennas), tau_.shape[0])
        amp_allan_dev = np.empty(shape=shape, dtype='f8')
        amp_allan_dev_error = np.empty(shape=shape, dtype='f8')
        phase_allan_dev = np.empty(shape=shape, dtype='f8')
        phase_allan_dev_error = np.empty(shape=shape, dtype='f8')
        for ia, a in enumerate(antennas):
            for i, t in enumerate(tau_):
                adev_, adev_err_, _ = adev(np.abs(gains[:, a]), dt, t)
                amp_allan_dev[ia, i] = adev_
                amp_allan_dev_error[ia, i] = adev_err_
                adev_, adev_err_, _ = adev(np.angle(gains[:, a]) *
                                           (180.0 / np.pi), dt, t)
                phase_allan_dev[ia, i] = adev_
                phase_allan_dev_error[ia, i] = adev_err_

        mean_amp_allen_dev = np.mean(amp_allan_dev, axis=0)
        mean_amp_allen_dev_err = np.mean(amp_allan_dev_error, axis=0)
        mean_phase_allen_dev = np.mean(phase_allan_dev, axis=0)
        mean_phase_allen_dev_err = np.mean(phase_allan_dev_error, axis=0)


        ax = fig.add_subplot(223)
        # for i, a in enumerate(antennas):
        #     ax.errorbar(tau_, amp_allan_dev[i, :], amp_allan_dev_error[i, :],
        #                 fmt='.-')
        ax.errorbar(tau_, mean_amp_allen_dev, mean_amp_allen_dev_err, fmt='.-')
        ax.set_xlabel('Tau [seconds]')
        ax.set_ylabel('Amplitude allan deviation')

        ax = fig.add_subplot(224)
        # for i, a in enumerate(antennas):
        #     ax.errorbar(tau_, phase_allan_dev[i, :],
        #                 phase_allan_dev_error[i, :], fmt='.-')
        ax.errorbar(tau_, mean_phase_allen_dev, mean_phase_allen_dev_err,
                    fmt='.-')
        ax.set_xlabel('Tau [seconds]')
        ax.set_ylabel('Phase allan deviation [degrees]')

        plt.show()


    # Verify in MATLAB with:
    #
    # >> load('fbm_gains.mat')
    # >> plot(1:800, abs(fbm_gains.gains(:, 1:10)))
    # >> plot(1:800, angle(fbm_gains.gains(:, 1:10)))
    #


# http://wiki.scipy.org/Cookbook/SignalSmooth
def smooth(x, window_len=11, window='hanning'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")
    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y



def adev(data, dt, tau):
    """Evaluate the Allan deviation of a time series.

    Args:
        data (np.array): Time series data.
        dt (float): Time series data increment, in seconds.
        tau (float): Averaging period, in seconds.

    Returns:
        The allan deviation of the series and error on the deviation.
    """
    rate = 1. / dt
    m = int(rate * tau)
    # Truncate to an even multiple of this tau value
    freq = data[0:len(data) - int(np.remainder(len(data), m))]
    f = np.reshape(freq, (m, -1), order='F')
    fa = np.mean(f, 0)
    fd = np.diff(fa)
    n = len(fa) - 1
    sm = np.sqrt(0.5 / n * (np.sum(fd**2)))
    sme = sm / np.sqrt(n)
    return sm, sme, n


def fbm(n, hurst):
    """Generate Fractional brownian motion noise.

    http://www.maths.uq.edu.au/~kroese/ps/MCSpatial.pdf

    Args:
        n (int): Length of the time series.
        H (float): Hurst parameter.

    Kwargs:
        seed (int): Random generator number seed.

    Returns:
        Time series array.
    """
    r = np.empty((n + 1,)) * np.NaN
    r[0] = 1
    for k in range(1, n + 1):
        a = 2.0 * hurst
        r[k] = 0.5 * ((k + 1)**a - 2 * k**a + (k - 1)**a)
    r = np.append(r, r[-2:0:-1])  # first row of circulant matrix
    lambda_ = np.real(np.fft.fft(r)) / (2 * n)  # Eigenvalues
    W = np.fft.fft(np.sqrt(lambda_) * (np.random.randn(2 * n) +
                                       1j * np.random.randn(2 * n)))
    W = n**(-hurst) * np.cumsum(np.real(W[0:n + 1]))  # Rescale
    return W


def eval_gain_amp(length, tstep, hurst, adev_fbm, sigma_wn, tau):
    """."""
    g_fbm = fbm(length, hurst)
    g_fbm *= adev_fbm / adev(g_fbm, tstep, tau)[0]
    g_wn = np.random.randn(length + 1,) * sigma_wn
    return (g_fbm + g_wn) + 1.0


def eval_gain_phase(n, tstep, H, adev_fbm, sigma_wn, tau):
    """."""
    p_fbm = fbm(n, H)
    p_fbm *= adev_fbm / adev(p_fbm, tstep, tau)[0]
    p_wn = np.random.randn(n + 1,) * sigma_wn
    phase = p_fbm + p_wn  # Combine FBM and WN signals
    return phase


def eval_complex_gain(n, dt, amp_H, amp_adev_fbm, amp_sigma_wn,
                      phase_H, phase_adev_fbm, phase_sigma_wn,
                      amp_std_t0, phase_std_t0, tau):
    """."""
    amp = eval_gain_amp(n, dt, amp_H, amp_adev_fbm, amp_sigma_wn, tau)
    phase = eval_gain_phase(n, dt, phase_H, phase_adev_fbm, phase_sigma_wn,
                            tau) * np.pi / 180.0
    amp_t0 = np.random.randn() * amp_std_t0
    phase_t0 = np.random.randn() * phase_std_t0 * (np.pi / 180.)
    gain = (amp + amp_t0) * np.exp(1.0j * (phase + phase_t0))
    return gain[0:n]


if __name__ == '__main__':
    main()


