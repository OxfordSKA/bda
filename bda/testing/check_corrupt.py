# -*- coding: utf-8 -*-
"""Quick test script for checking corrupted ms against each other.

Run with casa
"""

import os
import numpy as np

# ms1 = os.path.join('out_default', 'vis', 'model.ms')
# ms2 = os.path.join('out_new', 'vis', 'model.ms')
#
# tb.open(ms1)
# ms1_data = tb.getcol('DATA')
# tb.close()
#
# tb.open(ms2)
# ms2_data = tb.getcol('DATA')
# tb.close()
#
# print '=' * 60
# print 'MODEL'
# print ms1_data[0, 0, 0]
# print ms2_data[0, 0, 0]
# print ms1_data[0, 0, 5]
# print ms2_data[0, 0, 5]
# print np.max(ms1_data - ms2_data)
# print '=' * 60


ms1 = os.path.join('out_default', 'vis', 'corrupted.ms')
ms2 = os.path.join('out_new', 'vis', 'corrupted.ms')

tb.open(ms1)
ms1_data = tb.getcol('DATA')
tb.close()

tb.open(ms2)
ms2_data = tb.getcol('DATA')
tb.close()

print '=' * 60
print ms1_data[0, 0, 0]
print ms2_data[0, 0, 0]
print ms1_data[0, 0, 5]
print ms2_data[0, 0, 5]
print np.max(ms1_data - ms2_data)
print '=' * 60
