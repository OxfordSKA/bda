# -*- coding: utf-8 -*-
"""Module containing utility bda for BDA simulation steps."""

import numpy
import math


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
    freq = data[0:len(data) - int(numpy.remainder(len(data), m))]
    f = numpy.reshape(freq, (m, -1), order='F')
    fa = numpy.mean(f, 0)
    fd = numpy.diff(fa)
    n = len(fa) - 1
    sm = numpy.sqrt(0.5 / n * (numpy.sum(fd**2)))
    sme = sm / numpy.sqrt(n)
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
    r = numpy.empty((n + 1,)) * numpy.NaN
    r[0] = 1
    for k in range(1, n + 1):
        a = 2.0 * hurst
        r[k] = 0.5 * ((k + 1)**a - 2 * k**a + (k - 1)**a)
    r = numpy.append(r, r[-2:0:-1])  # first row of circulant matrix
    lambda_ = numpy.real(numpy.fft.fft(r)) / (2 * n)  # Eigenvalues
    W = numpy.fft.fft(numpy.sqrt(lambda_) * (numpy.random.randn(2 * n) +
                                             1j * numpy.random.randn(2 * n)))
    W = n**(-hurst) * numpy.cumsum(numpy.real(W[0:n + 1]))  # Rescale
    return W


def eval_gain_amp(length, tstep, hurst, adev_fbm, sigma_wn, tau):
    """."""
    g_fbm = fbm(length, hurst)
    g_fbm *= adev_fbm / adev(g_fbm, tstep, tau)[0]
    g_wn = numpy.random.randn(length + 1) * sigma_wn
    return (g_fbm + g_wn) + 1.0


def eval_gain_phase(n, tstep, H, adev_fbm, sigma_wn, tau):
    """."""
    p_fbm = fbm(n, H)
    p_fbm *= adev_fbm / adev(p_fbm, tstep, tau)[0]
    p_wn = numpy.random.randn(n + 1) * sigma_wn
    phase = p_fbm + p_wn  # Combine FBM and WN signals
    return phase


def eval_complex_gain(n, dt, amp_H, amp_adev_fbm, amp_sigma_wn,
                      phase_H, phase_adev_fbm, phase_sigma_wn,
                      amp_std_t0, phase_std_t0, tau):
    amp = eval_gain_amp(n, dt, amp_H, amp_adev_fbm, amp_sigma_wn, tau)
    phase = eval_gain_phase(n, dt, phase_H, phase_adev_fbm, phase_sigma_wn,
                            tau)
    phase = numpy.radians(phase)
    amp_t0 = numpy.random.randn() * amp_std_t0
    phase_t0 = numpy.random.randn() * math.radians(phase_std_t0)
    gain = (amp + amp_t0) * numpy.exp(1.0j * (phase + phase_t0))
    return gain[0:n]


def byteify(input):
    """Convert unicode string."""
    if isinstance(input, dict):
        return {byteify(key):byteify(value) for key,value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input
