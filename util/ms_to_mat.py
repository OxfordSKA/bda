# Run using:
# casapy --nologger --nogui --log2term -c ms_to_mat <MS name>

import os
import sys
import scipy.io
import collections

# Get MS name from command line.
ms_name = sys.argv[-1]

# Create top-level dictionary.
d = collections.OrderedDict()
k = os.path.splitext(os.path.basename(ms_name))[0]
d[k] = collections.OrderedDict()

# Get data from main table.
tb.open(ms_name)
col_names = tb.colnames()
sub_tables = tb.getkeywords()
for c in col_names:
    try:
        d[k][c] = tb.getcol(c)
    except:
        d[k][c] = []
tb.close()

# Get data from sub-tables.
for s in sorted(sub_tables):
    if str(sub_tables[s]).startswith('Table'):
        tb.open(ms_name + '/' + s)
        d[k][s] = collections.OrderedDict()
        col_names = tb.colnames()
        for c in col_names:
            try:
                d[k][s][c] = tb.getcol(c)
            except:
                d[k][s][c] = []
        tb.close()

# Write data to MATLAB file.
scipy.io.savemat(ms_name + '.mat', d, do_compression=True)
