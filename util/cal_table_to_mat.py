#!/usr/bin/python

import numpy as np
import scipy.io
import os

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
    # ----------------------------------------------------------
    mat_file = os.path.join('gains.mat')
    cal_table = os.path.join('vis', 'test.cal')
    # ----------------------------------------------------------

    gains, num_antennas, num_times, dt = load_gains(cal_table)
    gains = gains[0, 0, :]

    data = {}
    data['gains'] = np.reshape(gains, (num_antennas, num_times))
    data['num_antennas'] = num_antennas
    data['num_times'] = num_times
    data['delta_t'] = dt
    scipy.io.savemat(mat_file, data, appendmat=False, do_compression=True,
                     oned_as='column')

if __name__ == "__main__":
    main()
