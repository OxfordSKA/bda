#!/usr/bin/python -u

"""
Simulate a Measurement Set for use with BDA tests.
"""

def create_sky_model(sky_file, ra, dec, I):
    """Create an OSKAR sky model."""
    import os
    if not os.path.dirname(sky_file):
        os.mkdir(os.path.dirname(sky_file))
    fh = open(sky_file, 'w')
    for ra_, dec_, I_ in zip(ra, dec, I):
        fh.write('%.14f, %.14f, %.3f' % (ra_, dec_, I_))
    fh.close()


def dict_to_ini(settings_dict, ini):
    """Convert a dictionary of settings to and OSKAR settings ini file."""
    import os
    ini_dir = os.path.dirname(ini)
    if not ini_dir == "" and not os.path.isdir(ini_dir):
        os.makedirs(ini_dir)
    for group in sorted(settings_dict.keys()):
        for key in sorted(settings_dict[group].keys()):
            key_ = group+key
            value_ = settings_dict[group][key]
            set(ini, key_, value_)


def run_interferometer(ini, verbose=True):
    """Run the OSKAR interferometer simulator."""
    from subprocess import call
    if verbose:
        call(['oskar_sim_interferometer', ini])
    else:
        call(['oskar_sim_interferometer', '-q', ini])


def create_settings(ini_file, sky, ms):
    """Create simulation settings file."""
    import os
    if not os.path.isdir(os.path.dirname(ms)):
        os.mkdir(os.path.dirname(ms))
    s = {}
    s['simulator/'] = {
        'max_sources_per_chunk': 1,
        'double_precision': 'true',
        'keep_log_file': 'false'
    }
    s['sky/'] = {
        'oskar_sky_model/file': sky
    }
    s['observation/'] = {
        'start_frequency_hz': 700e6,
        'start_time_utc': 57086.113194,
        'length': 320.0,
        'num_time_steps': 200,
        'phase_centre_ra_deg': -90.3545848760,
        'phase_centre_dec_deg': -8.5711239906
    }
    s['telescope/'] = {
        'longitude_deg': 21.442909,
        'latitude_deg': -30.739475,
        'input_directory': 'models/SKA1_mid_combined_rmax_5.000_km.tm',
        'pol_mode': 'Scalar',
        'station_type': 'Isotropic beam'
    }
    s['interferometer/'] = {
        'time_average_sec': 1.6,
        'channel_bandwidth_hz': 4000,
        'ms_filename': ms
    }
    dict_to_ini(s, ini_file)
    return s

def oskar_sim():
    import os
    # ---------------------------------------------
    ini = os.path.join('ini', 'test.ini')
    ms = os.path.join('vis', 'test.ms')
    sky = os.path.join('models', 'sky.osm')
    ra0 = -90.3545848760
    dec0 = -8.5711239906
    # ---------------------------------------------
    create_sky_model(sky, [ra0], [dec0+0.9], [1.0])
    create_settings(ini, sky, ms)
    run_interferometer(ini, verbose=False)
    return ms

def main():
    import time
    t0 = time.time()
    print '+ MS simulation...'
    ms = oskar_sim()
    print '  - Finished simulation in %.3fs' % (time.time()-t0)
    print '  - MS : %s' % ms

if __name__ == "__main__":
    main()
