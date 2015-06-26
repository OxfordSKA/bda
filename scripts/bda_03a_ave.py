# -*- coding: utf-8 -*-
"""BDA with custom BDA code."""

import numpy as np
import os
import shutil
import time
import sys
import subprocess


if __name__ == "__main__":
    # -------------------------------------------------------------------------
    dt = 0.08  # Correlator dump time. TODO(BM) get this from the MS.
    idt_max = 100
    dt_max = '%.2fs' % (idt_max * dt)  # Maximum allowed averaging time.
    max_fact = 1.01   # Maximum amplitude loss factor.
    fov_radius = 0.9  # Field of view radius (input into mstransform)
    # -------------------------------------------------------------------------

    cmd = 'src/bda'

    ms = os.path.join('out_new', 'vis', 'model.ms')
    subprocess.call([cmd, ms,
                     '%.2f' % max_fact,
                     '%.1f' % fov_radius,
                     '%i' % idt_max])

    # ms = os.path.join('out_new', 'vis', 'corrupted.ms')
    # subprocess.call([cmd, ms,
    #                  '%.2f' % max_fact,
    #                  '%.1f' % fov_radius,
    #                  '%i' % idt_max])
