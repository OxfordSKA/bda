# -*- coding: utf-8 -*-
"""."""

import numpy as np
import sys
import os


def angular_distance(ra0, dec0, RA, Dec):
    """."""
    deltaRA = (RA - ra0) / 2.0
    deltaDec = (Dec - dec0) / 2.0
    sinDeltaRA = np.sin(deltaRA * (np.pi / 180.0))
    sinDeltaDec = np.sin(deltaDec * (np.pi / 180.0))
    sinDeltaDecSq = np.power(sinDeltaDec, 2.0)
    sinDeltaRASq = np.power(sinDeltaRA, 2.0)
    cosDec = np.cos(Dec * (np.pi / 180.0))
    cosDec0 = np.cos(dec0 * (np.pi / 180.0))
    arg = sinDeltaDecSq + (cosDec0 * cosDec) * sinDeltaRASq
    arg = np.sqrt(arg)
    r = 2.0 * np.arcsin(arg)
    return r * (180.0 / np.pi)


def create_region_file(filename, reg_file):
    """."""
    # Read the OSM file
    osm_data = np.loadtxt(filename)

    coords = osm_data[:, 0:2]

    if (os.path.exists(reg_file)):
        os.remove(reg_file)

    file = open(reg_file, 'w')
    file.write('global color=green dashlist=8 3 width=1 '
               'font="helvetica 10 normal roman" select=1 highlite=1 dash=0 '
               'fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
    for i in range(0, len(coords)):
        file.write('point(%10.6f, %10.6f) # point=x\n' %
                   (coords[i, 0], coords[i, 1]))
    file.close()


def create_region_file_with_filter(filename, ra0_deg, dec0_deg,
                                   radius_deg, reg_file):
    """."""
    # Read the OSM file
    osm_data = np.loadtxt(filename, delimiter=',')

    # Evaluate radius values
    r = angular_distance(ra0_deg, dec0_deg, osm_data[:, 0], osm_data[:, 1])
    data = np.zeros((osm_data.shape[0], osm_data.shape[1] + 1))
    data[:, 0] = r
    data[:, 1:] = osm_data
    keep = data[:, 0] <= radius_deg
    data = data[keep, :]
    coords = data[:, 1:3]

    if (os.path.exists(reg_file)):
        os.remove(reg_file)

    file = open(reg_file, 'w')
    file.write('global color=green dashlist=8 3 width=1 '
               'font="helvetica 10 normal roman" select=1 highlite=1 dash=0 '
               'fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
    for i in range(0, len(coords)):
        file.write('point(%10.6f, %10.6f) # point=x\n' %
                   (coords[i, 0], coords[i, 1]))
    file.close()


if __name__ == '__main__':

    if (not (len(sys.argv) == 2 or len(sys.argv) == 5)):
        print 'Usage: '
        print ('  osm_to_ds9_regions <osm file> <RA0[deg]> '
               '<Dec0[deg]> [radius[deg]]')
        exit(1)

    osm_file = sys.argv[1]

    if (len(sys.argv) > 2):
        ra0_deg = float(sys.argv[2])
        dec0_deg = float(sys.argv[3])
        radius_deg = float(sys.argv[4])

    reg_file = osm_file+'.reg'

    print '-' * 60
    print 'Input osm file =', osm_file
    if len(sys.argv) == 5:
        print 'RA0 [deg] =', ra0_deg
        print 'Dec0 [deg] =', dec0_deg
        print 'Radius [deg] =', radius_deg
    print 'Output region file =', reg_file
    print '-' * 60, '\n'

    if len(sys.argv) == 5:
        create_region_file_with_filter(osm_file, ra0_deg, dec0_deg,
                                       radius_deg, reg_file)
    else:
        create_region_file(osm_file, reg_file)
