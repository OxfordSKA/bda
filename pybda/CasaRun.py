# -*- coding: utf-8 -*-

import pexpect
import sys
import os

class CasaRun(object):

    def __init__(self):
        cmd = [ '--nologger',
                '--nogui',
                '--log2term',
                '--nologfile',
                '--colors=NoColor']
        self.child = None
        try:
            self.child = pexpect.spawn('casapy',
                                       cmd,
                                       cwd=os.path.abspath(os.path.curdir),
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
        self.child.sendline(cmd)
        self.child.expect(self.prompt)

    def run_script(self, script, vars={}):
        for name in vars:
            cmd = "{} = {}".format(name, vars[name])
            self.run_command(cmd)
        cmd = "execfile('{}')".format(script)
        self.run_command(cmd)

    def open_table(self, table):
        pass
