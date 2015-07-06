
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
    idt_max = 50
    max_fact = 1.002   # Maximum amplitude loss factor.
    fov_radius = 0.9  # Field of view radius for max_fact
    # -------------------------------------------------------------------------

    cmd = 'src/bda'
    ms = os.path.join(sim_dir, 'vis', 'model.ms')
    subprocess.call([cmd, ms,
                     '%.5f' % max_fact,
                     '%.3f' % fov_radius,
                     '%i' % idt_max])

    print ''
    print '*' * 60
    print '*' * 60
    print '*' * 60
    print ''

    cmd = 'src/bda_2'
    ms = os.path.join(sim_dir, 'vis', 'corrupted.ms')
    subprocess.call([cmd, ms,
                     '%.5f' % max_fact,
                     '%.3f' % fov_radius,
                     '%i' % idt_max])
