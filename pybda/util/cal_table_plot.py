#!/usr/bin/python

import os
import matplotlib.pyplot as plt
import numpy as np


def load_gains(cal_table):
    tb.open(cal_table, nomodify=True)
    gains = tb.getcol('CPARAM')
    tb.close()
    tb.open(cal_table+'/ANTENNA')
    num_antennas = tb.nrows()
    tb.close()
    num_times = gains.shape[2]/num_antennas
    tb.open(cal_table+'/OBSERVATION')
    time_range = tb.getcol('TIME_RANGE')
    tb.close()
    dt = (time_range[1][0]-time_range[0][0]) / num_times
    return gains, num_antennas, num_times, dt


def main():
    # ------------------------------------------------------------
    # cal_table = os.path.join('vis', 'corrupted.gains')
    cal_table = os.path.join('bench_02', 'calibrated.gains')
    plot_num_stations = 10  # Number of stations to plot gains for
    # ------------------------------------------------------------
    gains, num_antennas, num_times, dt = load_gains(cal_table)
    if plot_num_stations == -1:
        plot_num_stations = num_antennas

    gains = gains[0, 0, :]
    if 'corrupted' in cal_table :
        gains = 1.0/gains
    x = np.arange(0, num_times)*dt
    fig, axes = plt.subplots(4, 1, sharex=True, sharey=False, figsize=(12,10))
    # for i in np.random.randint(0, num_antennas, plot_num_stations):
    for i in range(0, plot_num_stations):
        axes[0].plot(x, np.abs(gains[i::num_antennas]))
        axes[1].plot(x, np.angle(gains[i::num_antennas])*(180.0/np.pi))
        axes[2].plot(x, np.real(gains[i::num_antennas]))
        axes[3].plot(x, np.imag(gains[i::num_antennas]))

    for axis in axes:
        axis.grid()

    # axes[0].set_title('%s : Gains for %i randomly selected stations' %
    #                   (cal_table, plot_num_stations))
    axes[0].set_title('%s : Gains for the first %i stations' %
                     (cal_table, plot_num_stations))
    axes[0].set_ylabel('gain amplitude')
    axes[1].set_ylabel('gain phase')
    axes[2].set_ylabel('real(gain)')
    axes[3].set_ylabel('imag(gain)')
    axes[3].set_xlabel('time [seconds]')
    axes[3].set_xlim(0, num_times*dt)

    plt.tight_layout()
    #plt.savefig(cal_table+'.png', transparent=True, frameon=False)
    plt.savefig(cal_table+'.png')
    plt.show()


if __name__ == "__main__":
    main()
