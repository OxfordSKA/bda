# -*- coding: utf-8 -*-
"""Module to evaluate corrupting gains for BDA simulations."""

from __future__ import (print_function, absolute_import)
import numpy
import time
import os


def allan_deviation(data, dt, tau):
    """
    Evaluate the Allan deviation of a time series.

    References:
        https://en.wikipedia.org/wiki/Allan_variance

    Args:
        data (array_like): Array of time series data.
        dt (float): Sample spacing of the time series data.
        tau (float): Interval at which to calculate the allan deviation.

    Returns:
        sm: Allan deviation
        sme: error on the allan deviation
        n: number of pairs in the Allan computation
    """
    data = numpy.asarray(data)
    num_points = data.shape[0]
    m = int(tau / dt)  # Number of samples in length tau
    data = data[:num_points - (num_points % m)]  # Resize to a multiple of m
    # Reshape into blocks of length m and take the average of each block.
    data = data.reshape((-1, m))
    data_mean = numpy.mean(data, axis=1)
    data_diff = numpy.diff(data_mean)
    n = data_diff.shape[0]
    a_dev = numpy.sqrt((0.5 / n) * (numpy.sum(data_diff**2)))
    a_dev_err = a_dev / numpy.sqrt(n)
    return a_dev, a_dev_err, n


def smooth(x, window_len=11, window='hanning'):
    x = numpy.asarray(x)
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")

    s = numpy.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='valid')
    return y[:x.shape[0]]


def fractional_brownian_motion(n, hurst):
    """Generate Fractional brownian motion noise.

    http://www.maths.uq.edu.au/~kroese/ps/MCSpatial.pdf

    Args:
        n (int): Length of the time series.
        hurst (float): Hurst parameter.

    Returns:
        Time series array.
    """
    a = 2.0 * hurst
    r = numpy.empty(n + 1) * numpy.NaN
    r[0] = 1.0
    k = numpy.arange(1, n + 1)
    r[1:] = 0.5 * ((k + 1)**a - 2.0 * k**a + (k - 1)**a)
    r = numpy.append(r, r[-2:0:-1])  # first row of circulant matrix
    l = numpy.real(numpy.fft.fft(r)) / (2 * n)
    w = numpy.fft.fft(numpy.sqrt(l) * (numpy.random.randn(2 * n) +
                                       1.0j * numpy.random.randn(2 * n)))
    w = n**(-hurst) * numpy.cumsum(numpy.real(w[:n+1]))  # Rescale
    return w[:n]


def eval_complex_gains(n, dt, hurst_amp, adev_amp, std_t_mid_amp,
                       hurst_phase, adev_phase, std_t_mid_phase,
                       smoothing_length=0, tau=1.0):
    amp = fractional_brownian_motion(n, hurst_amp)
    if smoothing_length > 0:
        amp = smooth(amp, smoothing_length)
    amp *= adev_amp / allan_deviation(amp, dt, tau)[0]
    amp += 1.0 - amp[n / 2] + (numpy.random.randn() * std_t_mid_amp)

    phase = fractional_brownian_motion(n, hurst_phase)
    if smoothing_length > 0:
        phase = smooth(phase, smoothing_length)
    phase *= adev_phase / allan_deviation(phase, dt, tau)[0]
    phase += -phase[n / 2] + (numpy.random.randn() * std_t_mid_phase)

    gain = amp * numpy.exp(1.0j * numpy.radians(phase))
    return gain


def allan_dev_spectrum(data, dt):
    data = numpy.asarray(data)
    tau_values = numpy.arange(2, data.shape[0] / 3, 5) * dt
    adev = numpy.empty_like(tau_values)
    adev_err = numpy.empty_like(tau_values)
    for i, tau in enumerate(tau_values):
        adev[i] = allan_deviation(data, dt, tau)[0]
        adev_err[i] = allan_deviation(data, dt, tau)[0]
    return adev, tau_values, adev_err


def test_unblocked():
    import matplotlib.pyplot as pyplot

    tau = 1.0
    num_steps = 6
    num_antennas = 50
    hurst_amp = 0.55
    adev_amp = numpy.linspace(1.e-5, 1.e-3, num_steps)
    std_t_mid_amp = 0.05

    hurst_phase = 0.55
    adev_phase = numpy.linspace(0.01, 0.2, num_steps)
    std_t_mid_phase = 5.0

    num_times = 5000
    dump_time = 0.1
    over_sample = 10
    smoothing_length = over_sample * 3

    n = num_times * over_sample
    dt = dump_time / float(over_sample)
    times = numpy.arange(n) * dt

    print('No. samples = %i' % n)
    print('No. steps = %i' % num_steps)

    gains = numpy.empty((num_steps, num_antennas, n), dtype='c16')

    t0 = time.time()
    for i in range(num_steps):
        print('%i %e' % (i, adev_amp[i]))
        for a in range(num_antennas):
            gains[i, a, :] = eval_complex_gains(n, dt, hurst_amp, adev_amp[i],
                                                std_t_mid_amp,
                                                hurst_phase, adev_phase[i],
                                                std_t_mid_phase,
                                                smoothing_length,
                                                tau)

    print('Time taken to generate gains = %.3f s' % (time.time() - t0))

    fig = pyplot.figure(figsize=(7, 7))
    fig.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95,
                        hspace=0.2, wspace=0.0)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    colors = ['r', 'g', 'b', 'm', 'y', 'k']
    print('plotting ...')
    y_max = 0.0
    for i in range(num_steps):
        for a in range(num_antennas):
            if i == num_steps - 1:
                ax1.plot(times, numpy.abs(gains[i, a, :]), '-', color=colors[i])
            s, t, _ = allan_dev_spectrum(numpy.abs(gains[i, a, :]), dt)
            ax2.semilogy(t, s, color=colors[i])
            ax2.set_xlim(0, tau * 3.0)
            y_max = max(y_max, s[numpy.argmax(t > tau * 3.0)])
            ax2.set_ylim(0, y_max * 1.05)

    y_max = numpy.max(gains[num_steps-1, :, :])
    print(numpy.abs(y_max))

    ax1.grid()
    ax1.set_xlabel('Observation length [seconds]', fontsize='small')
    ax1.set_ylabel('Gain amplitude', fontsize='small')
    ax2.grid()
    ax2.set_ylabel('Allan deviation', fontsize='small')
    ax2.set_xlabel('Sample period, tau [seconds]', fontsize='small')
    if os.path.isfile('gain_amp.png'):
        os.remove('gain_amp.png')
    pyplot.savefig('gain_amp.png')

    # =========================================================================
    fig = pyplot.figure(figsize=(7, 7))
    fig.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95,
                        hspace=0.2, wspace=0.0)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    colors = ['r', 'g', 'b', 'm', 'y', 'k']
    print('plotting ...')
    y_max = 0.0
    for i in range(num_steps):
        for a in range(num_antennas):
            if i == num_steps - 1:
                ax1.plot(times, numpy.degrees(numpy.angle(gains[i, a, :])),
                         '-', color=colors[i])
            s, t, _ = allan_dev_spectrum(
                numpy.degrees(numpy.angle(gains[i, a, :])), dt)
            ax2.semilogy(t, s, color=colors[i])
            ax2.set_xlim(0, tau * 3.0)
            y_max = max(y_max, s[numpy.argmax(t > tau * 3.0)])
            ax2.set_ylim(0, y_max * 1.05)

    y_max = numpy.max(gains[num_steps-1, :, :])
    print(numpy.abs(y_max))

    ax1.grid()
    ax1.set_xlabel('Observation length [seconds]', fontsize='small')
    ax1.set_ylabel('Gain phase', fontsize='small')
    ax2.grid()
    ax2.set_ylabel('Allan deviation', fontsize='small')
    ax2.set_xlabel('Sample period, tau [seconds]', fontsize='small')
    if os.path.isfile('gain_phase.png'):
        os.remove('gain_phase.png')
    pyplot.savefig('gain_phase.png')


if __name__ == '__main__':
    test_unblocked()

