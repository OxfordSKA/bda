# -*- coding: utf-8 -*-

import drivecasa
import os
import subprocess
import pexpect
import sys

class Casapy(object):

    def __init__(self):
        cmd = [ '--nologger',
                '--nogui',
                '--log2term',
                '--colors=NoColor']
        self.child = None
        try:
            self.child = pexpect.spawn('casapy',
                                       cmd,
                                       cwd=os.path.curdir,
                                       timeout=600)
            self.child.logfile_read = sys.stdout
            self.prompt = r'CASA <[0-9]+>:'
            self.child.expect(self.prompt, timeout=60)
        except pexpect.TIMEOUT:
            print 'failed to start CASAPY!!!!'
            self.child = None

        if self.child is None:
            raise RuntimeError('failed to create child :(')

    def run_command(self, cmd):
        # exec_cmd = "execfile('{}')".format(os.path.abspath(cmd))
        self.child.sendline(cmd)
        self.child.expect(self.prompt)

# casa = drivecasa.Casapy(working_dir=os.path.curdir,
#                         casa_logfile=False,
#                         commands_logfile='commands.log',
#                         echo_to_stdout=True)
# casa.run_script_from_file('test/casa_script.py')

#casa = drivecasa.Casapy()

casa = Casapy()

casa.run_command("print 'hello'")
#casa.run_command("execfile('test/casa_script.py', {'hello':2})")
casa.run_command("execfile('test/casa_script.py', globals(), {'ms_name': 'TEMP/model.ms'})")
