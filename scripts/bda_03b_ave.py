#!/usr/bin/python
"""BDA averaging without mstransform."""

import numpy as np
import os
import shutil
import time
from progressbar import ProgressBar, Percentage, Bar, ETA


def inv_sinc(arg):
    """Newton-Raphson method for calculating arcsinc(x), from Obit."""
    import numpy as np
    x1 = 0.001
    for i in range(0, 1000):
        x0 = x1
        a = x0 * np.pi
        x1 = x0 - ((np.sin(a) / a) - arg) / \
            ((a * np.cos(a) - np.pi * np.sin(a)) / (a**2))
        if (np.fabs(x1 - x0) < 1.0e-6):
            break
    return x1


def load_ms(ms):
    """Load a measurment set."""
    tb.open(ms)
    dtype = [('time', 'f8'), ('uvw', 'f8'), ('uu', 'f8'), ('vv', 'f8'),
             ('ww', 'f8'), ('a1', 'i4'), ('a2', 'i4'), ('interval', 'f8'),
             ('data', 'c16')]
    if 'MODEL_DATA' in tb.colnames():
        dtype.append(('model', 'c16'))
    if 'CORRECTED_DATA' in tb.colnames():
        dtype.append(('corrected', 'c16'))
    vis = np.zeros((tb.nrows(),), dtype=dtype)
    vis['time'] = tb.getcol('TIME')
    uvw_ = tb.getcol('UVW')
    vis['uu'] = uvw_[0, :]
    vis['vv'] = uvw_[1, :]
    vis['ww'] = uvw_[2, :]
    vis['a1'] = tb.getcol('ANTENNA1')
    vis['a2'] = tb.getcol('ANTENNA2')
    vis['interval'] = tb.getcol('INTERVAL')
    vis['data'] = tb.getcol('DATA')
    if 'MODEL_DATA' in tb.colnames():
        vis['model'] = tb.getcol('MODEL_DATA')
    if 'CORRECTED_DATA' in tb.colnames():
        vis['corrected'] = tb.getcol('CORRECTED_DATA')
    tb.close()
    return vis


def get_time_info(ms):
    """Extract time information from a MS."""
    tb.open(ms, nomodify=True)
    times = np.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [np.min(times), np.max(times)]
    num_times = len(times)
    length = time_range[1] - time_range[0]
    dt = length / (num_times - 1)
    return num_times, time_range, length, dt


def get_num_antennas(ms):
    """Return number of antennas in a MS."""
    tb.open(ms + '/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations


def get_delta_uv_max(ms, fov_radius, max_fact):
    """Evaluate maximum allowed baseline coordinate momvement.

    Evalaute the maximum uvw distance allowed to achieve a the specified
    amplitude drop for the specified radius.
    """
    tb.open(ms + '/SPECTRAL_WINDOW', nomodify=True)
    freqs = tb.getcell('CHAN_FREQ', 0)
    tb.close()
    if len(freqs) != 1:
        print 'ERROR: can only use single channel data.'
        return
    freq = freqs[0]
    wavelength = 299792458.0 / freq
    delta_uv = inv_sinc(1.0 / max_fact) / (fov_radius * (np.pi / 180.))
    delta_uv *= wavelength
    return delta_uv


def copy_ms(ms_in, ms_out):
    """Make a copy of a MS without copying rows of the main table."""
    tb.open(ms_in, nomodify=True)
    tbcopy = tb.copy(ms_out, deep=True, valuecopy=True, norows=True,
                     returnobject=False)
    if tbcopy:
        tbcopy.close()
    tb.close()
    tb.open(ms_in + '/ANTENNA')
    tbcopy = tb.copy(ms_out + '/ANTENNA', deep=True, valuecopy=True,
                     norows=False, returnobject=False)
    if tbcopy:
        tbcopy.close()
    tb.close()
    # TODO(BM) copy other subtables...


if __name__ == "__main__":
    # -------------------------------------------------------------------------
    dt = 1.6  # Correlator dump time. TODO(BM) get this from the MS.
    idt_max = 5
    dt_max = float(idt_max) * dt
    fov_radius = 0.9  # Field of view radius (input into mstransform)
    max_fact = 1.01  # Maximum amplitude factor by which a source can drop.
    # -------------------------------------------------------------------------

    t0 = time.time()
    ms_in = os.path.join('vis', 'model.ms')
    ms_out = os.path.join('vis', 'model_ave.ms')

    if os.path.isdir(ms_out):
        print '+ INFO: Removing existing averaged MS (%s).' % (ms_out)
        shutil.rmtree(ms_out)

    if os.path.isdir(ms_in):
        # =====================================================================
        # =====================================================================
        duv_max = get_delta_uv_max(ms_in, fov_radius, max_fact)
        copy_ms(ms_in, ms_out)

        # Load un-averaged visibility data into memory
        vis = load_ms(ms_in)
        num_antennas = get_num_antennas(ms_in)
        num_baselines = (num_antennas * (num_antennas - 1)) / 2
        num_times, time_range, length, dt = get_time_info(ms_in)
        print '+ Before averaging: %i rows == %i baselines x %i times' % \
            (len(vis), num_baselines, num_times)

        # per baseline temporary arrays.
        # TODO(BM) instead of writing one row at a time, write every time a
        # buffer of rows is filled (with variable buffer length).
        b_uu0 = np.zeros((num_baselines,), dtype=np.double)
        b_vv0 = np.zeros((num_baselines,), dtype=np.double)
        b_ww0 = np.zeros((num_baselines,), dtype=np.double)
        b_t0 = np.zeros((num_baselines,), dtype=np.double)
        ave_count = np.zeros((num_baselines,), dtype='i4')
        ave_data = np.zeros((num_baselines,), dtype='c16')
        ave_uvw = np.zeros((num_baselines,), dtype=('f8', 3))
        ave_time = np.zeros((num_baselines,), dtype='f8')

        for i in range(0, num_baselines):
            b_uu0[i] = vis['uu'][i]
            b_vv0[i] = vis['vv'][i]
            b_ww0[i] = vis['ww'][i]
            b_t0[i] = vis['time'][i]

        ant_pq = np.zeros((num_baselines,), dtype=('i4', 2))
        i = 0
        for p in range(0, num_antennas):
            for q in range(p + 1, num_antennas):
                ant_pq[i] = np.array((p, q))
                i += 1

        print '-' * 80
        print 'max fact : %f' % max_fact
        print 'duv_max  : %.3f' % duv_max
        print 'idt_max  : %i' % idt_max
        print 'dt_max   : %.3f' % dt_max
        print '-' * 80

        show_progress = True
        tb.open(ms_out, nomodify=False)

        if show_progress:
            widgets = [Percentage(), ' ', Bar(), ' ~', ETA()]
            pbar = ProgressBar(widgets=widgets, maxval=len(vis))
            pbar.start()

        for itime in range(0, num_times):
            for ib in range(0, num_baselines):

                idx = itime * num_baselines + ib

                if show_progress:
                    pbar.update(idx + 1)

                # Add to average.
                ave_data[ib] += vis['data'][idx]
                ave_uvw[ib][0] += vis['uu'][idx]
                ave_uvw[ib][1] += vis['vv'][idx]
                ave_uvw[ib][2] += vis['ww'][idx]
                ave_time[ib] += vis['time'][idx]
                ave_count[ib] += 1

                # Work out distances for next sample and decide if we need to
                # start a new average. The look ahead avoids having to
                # worry about edge effects ?
                if itime < num_times - 1:
                    idx1 = (itime + 1) * num_baselines + ib
                    b_duu = vis['uu'][idx1] - b_uu0[ib]
                    b_dvv = vis['vv'][idx1] - b_vv0[ib]
                    b_dww = vis['ww'][idx1] - b_ww0[ib]
                    b_dt = vis['time'][idx1] - b_t0[ib]
                    b_duvw = (b_duu**2 + b_dvv**2 + b_dww**2)**0.5
                    b_idt = np.round(b_dt / dt)

                if b_duvw >= duv_max or b_idt >= idt_max or \
                   itime == num_times - 1:

                    # write the row into the averaged ms
                    tb.addrows(1)
                    row = tb.nrows() - 1
                    # FIXME(BM) putting cells one at a time like this is
                    # just as expensive as mstransform, need to buffer in
                    # blocks
                    tb.putcell('UVW', row, ave_uvw[ib] / ave_count[ib])
                    tb.putcell('ANTENNA1', row, ant_pq[ib][0])
                    tb.putcell('ANTENNA2', row, ant_pq[ib][1])
                    tb.putcell('DATA', row, [[ave_data[ib] / ave_count[ib]]])
                    tb.putcell('EXPOSURE', row, ave_count[ib] * dt)
                    tb.putcell('INTERVAL', row, ave_count[ib] * dt)
                    t_ = ave_time[ib] / ave_count[ib]
                    tb.putcell('TIME', row, t_)
                    tb.putcell('TIME_CENTROID', row, t_)

                    # Append to the average data buffer and update inital time
                    # and baseline positions if not the last time.
                    if itime < num_times - 1:
                        idx0 = (itime + 1) * num_baselines + ib
                        b_uu0[ib] = vis['uu'][idx0]
                        b_vv0[ib] = vis['vv'][idx0]
                        b_ww0[ib] = vis['ww'][idx0]
                        b_t0[ib] = vis['time'][idx0]
                        ave_data[ib] = 0.0 + 0.0j
                        ave_uvw[ib][0] = 0.0
                        ave_uvw[ib][1] = 0.0
                        ave_uvw[ib][2] = 0.0
                        ave_time[ib] = 0.0
                        ave_count[ib] = 0

        tb.flush()
        if show_progress:
            pbar.finish()
        print 'Final number of rows = %i' % (tb.nrows())
        tb.close()

    print ''
    print '+ Time taken in averaging = %.3fs [%s]' % (time.time() - t0, ms_out)
