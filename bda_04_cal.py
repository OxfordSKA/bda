#!/usr/bin/python

import os
import time
import numpy as np


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
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

    # # -----------------------------------------
    # ms_in = os.path.join('vis', 'test_cor_ave.ms')
    # ms = os.path.join('vis', 'test_cor_ave_cal.ms')
    # cal_table = ms + '.gcal'
    # # -----------------------------------------

    # copytree(ms_in, ms)

    # t0 = time.time()
    # run_gaincal(ms, cal_table)
    # print '+ Gaincal completed in %.3fs' % (time.time()-t0)
    #
    # t0 = time.time()
    # run_applycal(ms, cal_table)
    # print '+ Applycal completed in %.3fs' % (time.time()-t0)
    #
    # print '+ Total time = %.3fs' % (time.time()-tAll)

    # -----------------------------------------
    ms_in = os.path.join('vis', 'test_cor.ms')
    ms = os.path.join('vis', 'test_cal.ms')
    cal_table = ms + '.gcal'
    # -----------------------------------------

    copytree(ms_in, ms)

    tAll = time.time()

    t0 = time.time()
    run_gaincal(ms, cal_table)
    print '+ Gaincal completed in %.3fs' % (time.time()-t0)

    t0 = time.time()
    run_applycal(ms, cal_table)
    print '+ Applycal completed in %.3fs' % (time.time()-t0)


if __name__ == "__main__":
    main()
