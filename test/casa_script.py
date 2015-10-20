# -*- coding: utf-8 -*-

def get_time_info(ms):
    """."""
    import numpy
    tb.open(ms, nomodify=True)
    times = numpy.unique(tb.getcol('TIME_CENTROID'))
    tb.close()
    time_range = [numpy.min(times), numpy.max(times)]
    num_times = len(times)
    length = time_range[1] - time_range[0]
    dt = length / (num_times - 1)
    return num_times, time_range, length, dt

# ms_name is passed as a local() via execfile().
print get_time_info(ms_name)
