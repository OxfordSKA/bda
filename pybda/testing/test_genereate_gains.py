# coding=utf-8
"""Script to test generation of complex gains."""

from __future__ import print_function, division, absolute_import
from pybda.utilities import (adev, fbm, eval_complex_gain,
                             eval_gain_amp, eval_gain_phase)
import matplotlib.pyplot as pyplot
import numpy
import time
from scipy import interpolate
import allantools


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
    s = numpy.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='valid')
    return y



def plot_gain_spectrum(complex_gains, dt):
    num_times = len(complex_gains[0])
    num_antennas = len(complex_gains)
    tau_ = numpy.arange(5, num_times / 3, 5) * dt
    print(num_antennas * tau_.shape[0])
    shape = (num_antennas, tau_.shape[0])
    amp_allan_dev = numpy.empty(shape=shape, dtype='f8')
    amp_allan_dev_error = numpy.empty(shape=shape, dtype='f8')
    phase_allan_dev = numpy.empty(shape=shape, dtype='f8')
    phase_allan_dev_error = numpy.empty(shape=shape, dtype='f8')
    t0 = time.time()
    for ia, a in enumerate(range(num_antennas)):
        for i, t in enumerate(tau_):
            amp = numpy.abs(complex_gains[a])
            # phase = numpy.degrees(numpy.angle(complex_gains[a]))
            phase = numpy.angle(complex_gains[a])
            adev_, adev_err_, _ = adev(amp, dt, t)
            amp_allan_dev[ia, i] = adev_
            amp_allan_dev_error[ia, i] = adev_err_

            adev_, adev_err_, _ = adev(phase, dt, t)
            phase_allan_dev[ia, i] = adev_
            phase_allan_dev_error[ia, i] = adev_err_
    print(' - Time taken = %.2f' % (time.time() - t0))

    mean_amp_allan_dev = numpy.mean(amp_allan_dev, axis=0)
    mean_phase_allan_dev = numpy.mean(phase_allan_dev, axis=0)

    # Errors using std.dev. of the Allan deviations.
    mean_amp_allan_dev_err = numpy.std(amp_allan_dev, axis=0)
    mean_phase_allan_dev_err = numpy.std(phase_allan_dev, axis=0)

    fig = pyplot.figure(figsize=(10.0, 8.0))
    fig.subplots_adjust(left=0.15, bottom=0.10, right=0.97, top=0.95,
                        hspace=0.1, wspace=0.0)
    x = numpy.arange(num_times) * dt
    ax = fig.add_subplot(111)
    ax.loglog(tau_, mean_amp_allan_dev, 'r-', label='amp')
    ax.loglog(tau_, mean_phase_allan_dev, 'b-', label='phase')
    ax.set_xlabel('Tau [seconds]', fontsize='small')
    ax.set_ylabel('Allan dev.', fontsize='small')
    ax.legend()
    pyplot.show()



def test1():
    """Generate gains for a number of antennas and evaluate their spectra.

    Check if the spectra can be manipulated in a predictable way by
    one of the parameters adev or sigma.
    """
    num_times = 2000
    dump_time = 0.1
    over_sample = 20
    n = num_times * over_sample  # Number of gain values (number of times)
    dt = dump_time / over_sample  # Time increment in data
    samples = 10
    niter = 5
    tau = 1.0
    hurst0 = 0.5
    allan_dev0 = 1.0e-4
    sigma = 0.0
    d_hurst = 0.0
    d_allan_dev = 1.0e-4

    fig = pyplot.figure(figsize=(18, 12))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    colors = ['r', 'g', 'b', 'k', 'y', 'm']
    times = numpy.arange(n) * dt
    tau_ = numpy.arange(2, num_times / 2, 2) * dt
    amp_adev = numpy.empty_like(tau_)

    for i in range(niter):
        hurst = hurst0 + d_hurst * i
        allan_dev = allan_dev0 + d_allan_dev * i
        print('%-3i H=%.2f adev=%.2e' % (i, hurst, allan_dev))
        for k in range(samples):
            gain_amp = eval_gain_amp(n, dt, hurst, allan_dev, sigma, tau)
            gain_amp = gain_amp[0:n]
            gains = smooth(gain_amp, window_len=over_sample)
            gain_amp = gains[0:n]
            for it, t in enumerate(tau_):
                amp_adev[it], _, _ = adev(gain_amp, dt, t)
                if i == 0 and k == 0 and it == 100:
                    print(t, adev(gain_amp, dt, t))
                    print(allantools.adev(gain_amp, 1./dt, t))
            ax1.plot(times, gain_amp, linestyle='-', marker='None',
                     markersize=5, color=colors[i])
            ax2.loglog(tau_, amp_adev, linestyle='-', marker='None',
                       markersize=5, color=colors[i])
            ax2.set_ylim(min(amp_adev.min(), ax2.get_ylim()[0]),
                         max(amp_adev.max(), ax2.get_ylim()[1]),)
        # ax1.grid()
        ax2.set_xlim(tau_.min(), tau_.max())
        # ax2.grid()
        ax2.plot(ax2.get_xlim(), [allan_dev, allan_dev], 'r--', lw=1)

    pyplot.show()


def test2():
    num_antennas = 20
    num_times = 300
    dump_time = 0.1
    over_sample = 10
    n = num_times * over_sample  # Number of gain values (number of times)
    dt = dump_time / over_sample  # Time increment in data
    hurst = 0.7
    adev = 1.0e-5
    sigma = 0.0
    tau = 1.0
    amp_std_t0 = 1.0
    phase_std_t0 = 0.0
    # FIXME-BM shift so midpoint (instead of start) has the value {}_std_t0
    t0 = time.time()
    gains = list()
    for i in range(num_antennas):
        gains.append(eval_complex_gain(n, dt, hurst, adev, sigma,
                                       hurst, adev, sigma, amp_std_t0,
                                       phase_std_t0, tau))
    print('- Time taken to generate gains = %.2fs' % (time.time() - t0))

    # gain = gains[0]
    # x = numpy.arange(n) * dt
    # fig = pyplot.figure()
    # ax = fig.add_subplot(211)
    # ax.plot(x, numpy.abs(gain), '-', marker='+', markersize=5)
    # ax.set_ylabel('amplitude')
    # ax.grid()
    # ax = fig.add_subplot(212)
    # ax.plot(x, numpy.degrees(numpy.angle(gain)), '-', marker='+', markersize=5)
    # ax.grid()
    # ax.set_ylabel('phase [degrees]')
    # ax.set_xlabel('time [seconds]')
    # pyplot.show()
    print('- Plotting spectrum ...')
    plot_gain_spectrum(gains, dt)


def test3():
    """Fit splines to the gain curves and evaluate their spectra"""


if __name__ == '__main__':
    test1()
