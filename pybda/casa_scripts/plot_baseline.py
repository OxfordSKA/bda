# coding=utf-8

from os.path import join
import numpy
import matplotlib.pyplot as plt

def _get_time_info(ms):
    """."""
    tb.open(ms, nomodify=True)
    times = numpy.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [numpy.min(times), numpy.max(times)]
    num_times = len(times)
    length = time_range[1] - time_range[0]
    dt = length / (num_times - 1)
    return num_times, time_range, length, dt


def _get_num_antennas(ms):
    tb.open(ms + '/ANTENNA', nomodify=True)
    num_stations = tb.nrows()
    tb.close()
    return num_stations

def _load_ms_data(ms):
    ms_in = join(sim_dir, 'calibrated.ms')
    tb.open(ms, nomodify=True)
    num_rows = tb.nrows()
    col_time = tb.getcol('TIME')
    col_uvw = tb.getcol('UVW')
    col_ant1 = tb.getcol('ANTENNA1')
    col_ant2 = tb.getcol('ANTENNA2')
    col_data = numpy.squeeze(tb.getcol('DATA'))
    col_model = numpy.squeeze(tb.getcol('MODEL_DATA'))
    col_corrected = numpy.squeeze(tb.getcol('CORRECTED_DATA'))
    dtype = [('time', 'f8'), ('a1', 'i8'), ('a2', 'i8'),
             ('uu', 'f8'), ('vv', 'f8'), ('ww', 'f8'),
             ('data', 'c16'), ('model', 'c16'), ('corrected', 'c16')]
    data = numpy.zeros((num_rows,), dtype=dtype)
    data['time'] = col_time
    data['a1'] = col_ant1
    data['a2'] = col_ant2
    data['uu'] = col_uvw[0, :]
    data['vv'] = col_uvw[0, :]
    data['ww'] = col_uvw[0, :]
    data['data'] = col_data
    data['model'] = col_model
    data['corrected'] = col_corrected
    tb.close()
    return data


def _plot(data, column, type, y_type, marker_color, marker, marker_edge_width,
          label=''):
    x = data['time']
    if y_type == 'amp':
        y = numpy.abs(data[column])
    elif y_type == 'phase':
        y = numpy.degrees(numpy.angle(data[column]))
    else:
        print 'Error, unknown y_type'
        return
    if marker[column] == 'o':
        marker_face_color = 'none'
    else:
        marker_face_color = marker_color[type]
    ax.plot(x, y,
            linestyle=':',
            color=marker_color[type],
            markeredgecolor=marker_color[type],
            markerfacecolor=marker_face_color,
            marker=marker[column],
            markeredgewidth=1.0,
            markersize=marker_size[column],
            label=label)
    # if y_type == 'amp':
    #     ax.set_ylim(1.00-0.0025, 1.00+0.0025)


def _plot_data(data, type, y_type, marker_color, marker, marker_edge_width,
               label=''):
    # _plot(data, 'data', type, y_type, marker_color, marker, marker_edge_width,
    #       'data ' + label)
    _plot(data, 'model', type, y_type, marker_color, marker, marker_edge_width,
          'model ' + label)
    # _plot(data, 'corrected', type, y_type, marker_color, marker, marker_edge_width,
    #       'corrected ' + label)


# =============================================================================
reference_data = _load_ms_data(join(sim_dir, 'calibrated.ms'))
expanded_bda_data = _load_ms_data(join(sim_dir, 'expanded_calibrated.ms'))
bda_data = _load_ms_data(join(sim_dir, 'bda_calibrated.ms'))

_, time_range, _, dt = _get_time_info(join(sim_dir, 'calibrated.ms'))
t0 = time_range[0] - dt/2
reference_data['time'] -= t0
expanded_bda_data['time'] -= t0
bda_data['time'] -= t0

marker_color = {
    'expanded': 'b',
    'ref': 'g',
    'bda': 'r'
}
marker = {
    'model': 'x',
    'corrected': 'o',
    'data': 's'
}
marker_size = {
    'model': 8.0,
    'corrected': 6.0,
    'data': 3.0
}
marker_edge_width = {
    'model': 1.5,
    'corrected': 1.5,
    'data': 1.0
}


for b in range(65, 67):
    baseline = (0, b)
    b_data = reference_data[reference_data['a1'] == baseline[0]]
    b_data = b_data[b_data['a2'] == baseline[1]]
    b_bda_data = bda_data[bda_data['a1'] == baseline[0]]
    b_bda_data = b_bda_data[b_bda_data['a2'] == baseline[1]]
    b_expanded_bda_data = expanded_bda_data[expanded_bda_data['a1'] == baseline[0]]
    b_expanded_bda_data = b_expanded_bda_data[b_expanded_bda_data['a2'] == baseline[1]]

    baseline_length = (b_data['uu']**2 + b_data['vv']**2 + b_data['ww']**2)**0.5
    baseline_length = numpy.mean(baseline_length)

    # =============================================================================
    # TODO-BM plot a delta-time axis by subtracting the correct midpoint
    # adjusted start time

    fig, axes2d = plt.subplots(nrows=2, sharex=True, figsize=(14, 10))
    fig.subplots_adjust(left=0.08, bottom=0.10, right=0.97, top=0.95,
                        hspace=0.18, wspace=0.0)

    # ==== Amplitude ========
    ax = axes2d[0]
    _plot_data(b_data, 'ref', 'amp', marker_color, marker,
               marker_edge_width)
    _plot_data(b_expanded_bda_data, 'expanded', 'amp', marker_color, marker,
               marker_edge_width)
    # _plot_data(b_bda_data, 'bda', 'amp', marker_color, marker,
    #            marker_edge_width)
    ax.set_title('baseline: %i - %i, length: %.2f m' % (baseline[0],
                                                        baseline[1],
                                                        baseline_length),
                 fontsize='small')
    ax.set_ylabel('Visibility amplitude [Jy]', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='small')
    ax.tick_params(axis='both', which='major', labelsize='small')

    # ==== Phase ========
    ax = axes2d[1]
    _plot_data(b_data, 'ref', 'phase', marker_color, marker,
               marker_edge_width, label='')
    _plot_data(b_expanded_bda_data, 'expanded', 'phase', marker_color, marker,
               marker_edge_width,
               label='expanded bda')
    # _plot_data(b_bda_data, 'bda', 'phase', marker_color, marker,
    #            marker_edge_width,
    #            label='bda')
    ax.set_xlabel('Time [seconds]', fontsize='small')
    ax.set_ylabel('Visibility phase [degrees]', fontsize='small')
    ax.grid(True)
    ax.tick_params(axis='both', which='minor', labelsize='small')
    ax.tick_params(axis='both', which='major', labelsize='small')
    ax = axes2d[1]
    ax.legend(bbox_to_anchor=(0, 1.02, 1.00, 0.5),
              loc=4,
              labelspacing=0.1,
              mode='expand',
              borderaxespad=0,
              handlelength=5.0,
              ncol=3,
              prop={'size':'small'})


    plt.savefig('vis_baseline_%03i_%03i.png' % (baseline[0], baseline[1]), dpi=300)
    #plt.show()
    # plt.close()
