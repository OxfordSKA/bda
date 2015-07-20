# -*- coding: utf-8 -*-
"""Simulate a Measurement Set for use with BDA tests."""

import os
import sys
import collections
import time
import subprocess
import numpy as np


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


def create_settings(ini_file, sky, telescope, ms, ra0, dec0):
    """Create simulation settings file."""
    if not os.path.isdir(os.path.dirname(ms)):
        os.mkdir(os.path.dirname(ms))
    # --------------------------------------------------------------
    # dt = 0.08  # seconds
    # num_times = 200
    dt = 0.02  # seconds
    num_times = 2000
    freq = 700.0e6  # Hz
    start_time = 57086.113194  # MJD UTC
    lon0 = 21.442909  # deg
    lat0 = -30.739475  # deg
    channel_bw = 0.0  # Hz
    # --------------------------------------------------------------
    s = collections.OrderedDict()
    s['simulator/'] = {
        # 'max_sources_per_chunk': 1,
        'double_precision': 'true',
        'keep_log_file': 'false'
    }
    s['sky/'] = {
        'oskar_sky_model/file': sky
    }
    s['observation/'] = {
        'start_frequency_hz': freq,
        'num_channels': 1,
        'start_time_utc': start_time,
        'length': num_times * dt,
        'num_time_steps': num_times,
        'phase_centre_ra_deg': ra0,
        'phase_centre_dec_deg': dec0
    }
    s['telescope/'] = {
        'longitude_deg': lon0,
        'latitude_deg': lat0,
        'input_directory': telescope,
        'pol_mode': 'Scalar',
        'station_type': 'Isotropic beam'
    }
    s['interferometer/'] = {
        'time_average_sec': dt,
        'channel_bandwidth_hz': channel_bw,
        'ms_filename': ms
    }
    dict_to_ini(s, ini_file)
    return s


def source_ring(ra0_deg, dec0_deg, radius_deg=0.9):
    """Generate a ring of sources around the phase centre."""
    # Positions of sources around the circle.
    pos = np.array([0, 70, 135, 135.127, 200, 270], dtype='f8') - 45
    pos *= np.pi / 180.0
    radius_lm = np.sin(radius_deg * np.pi / 180.0)
    l = radius_lm * np.cos(pos)
    m = radius_lm * np.sin(pos)

    # Project onto sphere at given position.
    dec0 = dec0_deg * np.pi / 180.0
    c = np.arcsin(np.sqrt(l**2 + m**2))
    dec = np.arcsin(np.cos(c) * np.sin(dec0) + m * np.cos(dec0))
    ra = np.arctan2(l, np.cos(dec0) * np.cos(c) - m * np.cos(dec0))
    return ra0_deg + ra * 180.0 / np.pi, dec * 180.0 / np.pi


def oskar_sim(sim_dir):
    """Run the OSKAR simulation."""
    # ---------------------------------------------
    ini = os.path.join(sim_dir, 'ini', 'test.ini')
    ms = os.path.join(sim_dir, 'vis', 'model.ms')
    sky = os.path.join(sim_dir, 'models', 'sky.osm')
    telescope = os.path.join('models', 'ska1_meerkat_mid_combined_july_2015.tm')
    ra0 = -90.3545848760  # deg
    dec0 = -8.5711239906  # deg
    # ---------------------------------------------
    ra, dec = source_ring(ra0, dec0)
    # create_sky_model(sky, [ra0], [dec0 + 0.9], [1.0])
    create_sky_model(sky, ra, dec, np.ones(ra.shape))
    create_settings(ini, sky, telescope, ms, ra0, dec0)
    run_interferometer(ini)
    return ms


if __name__ == "__main__":
    if len(sys.argv) - 1 != 1:
        print 'Usage:'
        print '  $ python scripts/bda_01_sim.py <simulation dir>'
        sys.exit(1)

    t0 = time.time()
    sim_dir = sys.argv[1]
    print '+ MS simulation...'
    ms = oskar_sim(sim_dir)
    print '  - Finished simulation in %.3fs' % (time.time() - t0)
    print '  - MS : %s' % ms
