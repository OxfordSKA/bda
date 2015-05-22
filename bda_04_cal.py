#!/usr/bin/python

import os
import time
import numpy as np


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            # shutil.copytree(s, d, symlinks, ignore)
            shutil.copytree(s, d)
        else:
            shutil.copy2(s, d)


def run_gaincal(ms, cal_table):
    gaincal(vis=ms, caltable=cal_table, field='', spw='', selectdata=False,
            solint='int', refant='0', minblperant=1, minsnr=3, gaintype='G',
            calmode='ap', append=False)


def run_applycal(ms, cal_table):
    # TODO(BM) replace with own function tested in test_apply_gains.py
    applycal(vis=ms, field='', spw='', selectdata=False,
            gaintable=[cal_table], gainfield=[''], interp=['nearest'],
            calwt=[False], applymode='calonly', flagbackup=False)


def main():
    tAll = time.time()

    
    # -------------------------------------------------------------------------
    ms_in = os.path.join('vis', 'corrupted.ms')
    ms = os.path.join('vis', 'calibrated.ms')
    cal_table = ms[:-3] + '.gains'
    # -------------------------------------------------------------------------

    if os.path.isdir(ms):
        print 'Removing existing MS : %s' % ms
        shutil.rmtree(ms)

    copytree(ms_in, ms)

    tAll = time.time()

    t0 = time.time()
    run_gaincal(ms, cal_table)
    print '+ Gaincal completed in %.3fs' % (time.time()-t0)

    t0 = time.time()
    run_applycal(ms, cal_table)
    print '+ Applycal completed in %.3fs' % (time.time()-t0)


    # -------------------------------------------------------------------------
    ms_in = os.path.join('vis', 'corrupted_ave.ms')
    ms = os.path.join('vis', 'calibrated_ave.ms')
    cal_table = ms[:-3] + '.gains'
    # -------------------------------------------------------------------------

    if os.path.isdir(ms):
        print 'Removing existing MS : %s' % ms
        shutil.rmtree(ms)

    print 'Copying %s to %s' % (ms_in, ms)
    os.system('cp -r %s %s' % (ms_in, ms)) 
    # copytree(ms_in, ms) # FIXME(BM) This command is Broken on SKA1 ... ?

    tAll = time.time()

    t0 = time.time()
    run_gaincal(ms, cal_table)
    print '+ Gaincal completed in %.3fs' % (time.time()-t0)

    t0 = time.time()
    run_applycal(ms, cal_table)
    print '+ Applycal completed in %.3fs' % (time.time()-t0)


if __name__ == "__main__":
    main()
    
    
