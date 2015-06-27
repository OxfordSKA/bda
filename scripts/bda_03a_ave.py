
"""BDA with custom BDA code."""

import os
import sys
import subprocess


if __name__ == "__main__":
    if len(sys.argv) - 1 < 1:
        print 'Usage:'
        print ('  $ python scripts/bda_03a_ave.py '
               '<simulation dir>')
        sys.exit(1)

    sim_dir = sys.argv[-1]
    if not os.path.isdir(sim_dir):
        print 'ERROR: simulation directory not found!'
        sys.exit(1)

    # -------------------------------------------------------------------------
    # dt = 0.1  # Correlator dump time. TODO(BM) get this from the MS.
    dt = 0.08  # Correlator dump time. TODO(BM) get this from the MS.
    idt_max = 100
    max_fact = 1.01   # Maximum amplitude loss factor.
    fov_radius = 0.9  # Field of view radius (input into mstransform)
    # -------------------------------------------------------------------------

    cmd = 'src/bda'
    ms = os.path.join(sim_dir, 'vis', 'model.ms')
    subprocess.call([cmd, ms,
                     '%.2f' % max_fact,
                     '%.1f' % fov_radius,
                     '%i' % idt_max])

    cmd = 'src/bda_2'
    ms = os.path.join(sim_dir, 'vis', 'corrupted.ms')
    subprocess.call([cmd, ms,
                     '%.2f' % max_fact,
                     '%.1f' % fov_radius,
                     '%i' % idt_max])
