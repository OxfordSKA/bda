# -*- coding: utf-8 -*-
"""Simulate a Measurement Set for use with BDA tests."""

import os
from os.path import join
from collections import OrderedDict
import subprocess
import numpy
import math


def create_sky_model(sky_file, ra, dec, stokes_i_flux):
    """Create an OSKAR sky model.

    Args:
        sky_file (string): filename path of the sky model being created.
        ra (float list): Right ascension, in degrees, of sources to put in the
            sky model.
        dec (float list): Declination, in degrees, of sources to put in the
            sky model.
        stokes_i_flux (float list): Stokes-I flux, in Jy, of sources to put in
            the sky model.
    """
    if not os.path.isdir(os.path.dirname(sky_file)):
        os.makedirs(os.path.dirname(sky_file))
    fh = open(sky_file, 'w')
    for ra_, dec_, I_ in zip(ra, dec, stokes_i_flux):
        fh.write('%.14f, %.14f, %.3f\n' % (ra_, dec_, I_))
    fh.close()


def source_ring(ra0_deg, dec0_deg, radius_deg=0.9):
    """Generate a ring of sources around the phase centre."""
    # Positions of sources around the circle.
    pos = numpy.array([0, 70, 135, 135.127, 200, 270, 135.-0.15], dtype='f8')
    pos -= 45.0
    pos *= math.pi / 180.0
    radius_lm = math.sin(radius_deg * math.pi / 180.0)
    l = radius_lm * numpy.cos(pos)
    m = radius_lm * numpy.sin(pos)

    # Project onto sphere at given position.
    dec0 = dec0_deg * math.pi / 180.0
    c = numpy.arcsin(numpy.sqrt(l**2 + m**2))
    dec = numpy.arcsin(numpy.cos(c) * math.sin(dec0) + m * math.cos(dec0))
    ra = numpy.arctan2(l, math.cos(dec0) * numpy.cos(c) - m * math.cos(dec0))
    return ra0_deg + ra * 180.0 / math.pi, dec * 180.0 / math.pi


def dict_to_ini(settings_dict, ini):
    """Convert a dictionary of settings to and OSKAR settings ini file."""
    ini_dir = os.path.dirname(ini)
    if not ini_dir == "" and not os.path.isdir(ini_dir):
        os.makedirs(ini_dir)
    for group in sorted(settings_dict):
        for key in sorted(settings_dict[group]):
            key_ = group + key
            value_ = settings_dict[group][key]
            subprocess.call(["oskar_settings_set", "-q", ini,
                            key_, str(value_)])


def run_interferometer(ini, verbose=True):
    """Run the OSKAR interferometer simulator."""
    if verbose:
        subprocess.call(["oskar_sim_interferometer", ini])
    else:
        subprocess.call(["oskar_sim_interferometer", "-q", ini])


def create_settings(settings):
    """Create simulation settings file."""
    sim = settings['sim']
    obs = sim['observation']
    tel = sim['telescope']

    ms_path = join(settings['path'], sim['output_ms'])
    if not os.path.isdir(os.path.dirname(ms_path)):
        os.mkdir(os.path.dirname(ms_path))

    sky_file = join(settings['path'], sim['sky_file'])

    s = OrderedDict()
    s['simulator/'] = {
        'double_precision': 'true',
        'keep_log_file': 'false'
    }
    s['sky/'] = {
        'oskar_sky_model/file': sky_file
    }
    s['observation/'] = {
        'start_frequency_hz': obs['freq_hz'],
        'num_channels': 1,
        'start_time_utc': obs['start_time_mjd'],
        'length': obs['num_times'] * obs['dt_s'],
        'num_time_steps': obs['num_times'],
        'phase_centre_ra_deg': obs['ra_deg'],
        'phase_centre_dec_deg': obs['dec_deg']
    }
    s['telescope/'] = {
        'longitude_deg': tel['lon_deg'],
        'latitude_deg': tel['lat_deg'],
        'input_directory': tel['path'],
        'pol_mode': 'Scalar',
        'station_type': 'Isotropic beam'
    }
    s['interferometer/'] = {
        'time_average_sec': obs['dt_s'],
        'channel_bandwidth_hz': obs['channel_bw_hz'],
        'ms_filename': ms_path
    }
    ini_file = join(settings['path'], sim['ini_file'])
    dict_to_ini(s, ini_file)
    return s


def run(settings, verbose=True):
    """Run the OSKAR simulation."""
    obs = settings['sim']['observation']
    ra, dec = source_ring(obs['ra_deg'], obs['dec_deg'])
    stokes_i = numpy.ones(ra.shape)
    stokes_i[-1] = 0.5
    # create_sky_model(sky, [ra0], [dec0 + 0.9], [1.0])
    sky_file = join(settings['path'], settings['sim']['sky_file'])
    create_sky_model(sky_file, ra, dec, stokes_i)
    create_settings(settings)
    ini_file = join(settings['path'], settings['sim']['ini_file'])
    run_interferometer(ini_file, verbose)
    return settings['sim']['output_ms']

