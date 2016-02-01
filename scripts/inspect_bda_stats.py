# -*- coding: utf-8 -*-
"""
Script to inspect a BDA measurement set to get a view of the distrubtion
of baselines that have been averaged.
"""

# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy


if __name__ == "__main__":
    # ms = 'TEST_gains_low_bda/calibrated_bda.ms'
    # ms = 'TEST_gains_high_bda/calibrated_bda.ms'
    ms = 'TEST_gains/calibrated_bda.ms'
    tb.open(ms, nomodify=True)
    col_uvw = tb.getcol('UVW')
    col_weight = numpy.squeeze(tb.getcol('WEIGHT'))
    tb.close()

    uu = col_uvw[0, :]
    vv = col_uvw[1, :]
    ww = col_uvw[1, :]
    r_3d = numpy.sqrt(numpy.sum(col_uvw**2, axis=0)) / 1.e3
    r_2d = numpy.sqrt(numpy.sum(col_uvw[0:2,:]**2, axis=0)) / 1.e3
    print r_2d.shape
    print col_weight.shape
    fig = plt.figure(figsize=(6.5, 5.0))
    # ax = fig.add_subplot(111, projection='3d')
    ax = fig.add_subplot(111)
    ax.plot(r_3d, col_weight, 'o')
    ax.set_xlabel('baseline length [km]')
    ax.set_ylabel('weight')
    ax.set_ylim(numpy.min(col_weight)-0.2, numpy.max(col_weight) + 0.2)
    # sc = ax.scatter(uu, vv, ww, c=col_weight,
    #                 marker='o', edgecolor='')
    # plt.colorbar(sc)
    plt.show()


