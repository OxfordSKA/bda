# -*- coding: utf-8 -*-
"""Module to evaluate maximum uv distance for averaging."""

import math
import argparse


def max_uv_distance(max_fact, radius_deg, freq_hz):
    """Function to evaluate maximum uv distance for given amplitude drop.

    Follows equation 1 from:
    ftp://ftp.cv.nrao.edu/NRAO-staff/bcotton/Obit/BLAverage.pdf

    Args:
        max_fact (float) : amplitude reduction factor
        radius_deg (float) : Source / FoV radius
        freq_hz (float) : Observation frequency, in Hz

    Returns:
        Maxium distance in the UV plane over which data can be averaged,
        in metres.
    """
    def inv_sinc(arg):
        """Newton-Raphson method for calculating arcsinc(x), from Obit."""
        x1 = 0.001
        for i in range(0, 1000):
            x0 = x1
            a = x0 * math.pi
            x1 = x0 - ((math.sin(a) / a) - arg) / \
                ((a * math.cos(a) - math.pi * math.sin(a)) / (a**2))
            if math.fabs(x1 - x0) < 1.0e-6:
                break
        return x1

    delta_uv = inv_sinc(1.0 / max_fact) / (radius_deg * (math.pi / 180.))
    wavelength = 299792458.0 / freq_hz
    delta_uv *= wavelength  # convert to metres
    return delta_uv

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Obtain maxium uv averaging '
                                                 'distance',
                                     epilog='''
                                     Example:
                                       $ python uv_distance.py 1.002 0.9 700e6
                                     ''')
    parser.add_argument('max_fact', help='Amplitude reduction factor., '
                                         'eg. 1.002', type=float)
    parser.add_argument('radius_deg', help='Distance at which max_fact '
                                           'applies. ie. source radius or '
                                           'edge of the FoV', type=float)
    parser.add_argument('freq_hz', help='Observation frequency, in Hz',
                        type=float)
    args = parser.parse_args()

    print '-' * 80
    print 'max_fact    : %.5f' % args.max_fact
    print 'radius, deg : %.3f' % args.radius_deg
    print 'freq, Hz    : %.3f' % args.freq_hz
    print '-' * 80

    max_uvw_distance = max_uv_distance(max_fact=args.max_fact,
                                       radius_deg=args.radius_deg,
                                       freq_hz=args.freq_hz)

    print 'Max uvw distance : %.6f m' % max_uvw_distance
