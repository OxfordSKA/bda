# Run using:
# casapy --nologger --nogui --log2term -c ms_to_mat <MS name>

import os
import sys
import scipy.io
import collections

# Get MS name from command line.
ms_path = os.path.abspath(sys.argv[-1])

if not os.path.isdir(ms_path):
    raise ValueError('Specified measurement set not found!.')

print '=' * 80
print 'Converting MS: %s' % ms_path
print '=' * 80

# Create top-level dictionary.
ms_name = os.path.splitext(os.path.basename(ms_path))[0]
ms_dict = collections.OrderedDict()

# Get data from main table.
tb.open(ms_path)
for col_name in tb.colnames():
    try:
        ms_dict[col_name] = tb.getcol(col_name)
    except RuntimeError:  # Thrown for empty columns
        ms_dict[col_name] = list()
# Get list of sub-tables.
sub_tables = tb.getkeywords()
tb.close()

# Get data from sub-tables.
for sub_table in sorted(sub_tables):
    if str(sub_tables[sub_table]).startswith('Table'):
        tb.open(os.path.join(ms_path, sub_table))
        ms_dict[sub_table] = collections.OrderedDict()
        for col_name in tb.colnames():
            try:
                ms_dict[sub_table][col_name] = tb.getcol(col_name)
            except RuntimeError:  # Thrown for empty columns
                ms_dict[sub_table][col_name] = list()
        tb.close()

# Write data to MATLAB file.
scipy.io.savemat(ms_path + '.mat', ms_dict, do_compression=True)
print '=' * 80
print 'MS conversion complete.'
print 'MAT file = %s' % ms_path + '.mat'
print '=' * 80
