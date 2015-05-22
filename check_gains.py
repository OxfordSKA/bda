#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as pp
import time

def load_gains(cal_table):
    dtype=[('time',  'f8'), ('a1', 'i4'), ('a2', 'i4'), ('cparam', 'c16', 2),
           ('paramerr', 'f8', 2), ('flag', 'i4',2), ('snr', 'f8', 2)]
    tb.open(cal_table, nomodify=True)
    gains = np.zeros((tb.nrows(),), dtype=dtype)
    gains['time'] = tb.getcol('TIME')
    gains['a1'] = tb.getcol('ANTENNA1')
    gains['a2'] = tb.getcol('ANTENNA2')
    gains_cparam = tb.getcol('CPARAM')
    gains['cparam'][:, 0] = gains_cparam[0, 0, :]
    gains['cparam'][:, 1] = gains_cparam[1, 0, :]
    gains_paramerr = tb.getcol('PARAMERR')
    gains['paramerr'][:, 0] = gains_paramerr[0, 0, :]
    gains['paramerr'][:, 1] = gains_paramerr[1, 0, :]
    gains_flag = tb.getcol('FLAG')
    gains['flag'][:, 0] = gains_flag[0, 0, :]
    gains['flag'][:, 1] = gains_flag[1, 0, :]
    gains_snr = tb.getcol('SNR')
    gains['snr'][:, 0] = gains_snr[0, 0, :]
    gains['snr'][:, 1] = gains_snr[1, 0, :]
    tb.close()

    tb.open(cal_table+'/ANTENNA')
    num_antennas = tb.nrows()
    tb.close()

    return gains, num_antennas

load_data = False
if load_data:
    # ----------------------------------------------------------------------
    cal_table_cal_ave = os.path.join('vis', 'calibrated_ave.gains')
    cal_table_cal = os.path.join('vis', 'calibrated.gains')
    cal_table_cor = os.path.join('vis', 'corrupted.gains')
    # ----------------------------------------------------------------------

    t0 = time.time()
    gains_cal_ave, num_antennas = load_gains(cal_table_cal_ave)
    print '+ Load of gain table complete in %.3f s.' % (time.time()-t0)
    print '+ No. antennas %i' % (num_antennas)

    t0 = time.time()
    gains_cal, num_antennas = load_gains(cal_table_cal)
    print '+ Load of gain table complete in %.3f s.' % (time.time()-t0)
    print '+ No. antennas %i' % (num_antennas)

    t0 = time.time()
    gains_cor, num_antennas = load_gains(cal_table_cor)
    print '+ Load of gain table complete in %.3f s.' % (time.time()-t0)
    print '+ No. antennas %i' % (num_antennas)

# Obtain gains for antenna a1 with phase referene of antenna a2
a1 = 100
a2 = 0
gains_cal_ave_a = gains_cal_ave[gains_cal_ave['a1']==a1]
gains_cal_ave_a = gains_cal_ave_a[gains_cal_ave_a['a2']==a2]
gains_cal_a = gains_cal[gains_cal['a1']==a1]
gains_cal_a = gains_cal_a[gains_cal_a['a2']==a2]
a2 = -1
gains_cor_a = gains_cor[gains_cor['a1']==a1]
gains_cor_a = gains_cor_a[gains_cor_a['a2']==a2]

# Remove flagged antenna
gains_cal_ave_a = gains_cal_ave_a[gains_cal_ave_a['flag'][:,0]==0]
gains_cal_a = gains_cal_a[gains_cal_a['flag'][:,0]==0]


gain_check = gains_cal_a['cparam'][:,0] * gains_cor_a['cparam'][:,0]
print 'gain amp check: amp %.3e+/-%.3e' % (np.mean(np.abs(gain_check))-1.0,
                                           np.std(np.abs(gain_check)))
print 'gain amp check: phase %.3e+/-%.3e deg' % (np.mean(np.angle(gain_check))*(180./np.pi),
                                                 np.std(np.angle(gain_check))*(180./np.pi))



# ave_ref = np.unique(gains_cal_ave['a2'])
# for a in ave_ref:
#     print a, np.sum(gains_cal_ave['a2']==a)/float(num_antennas)
#
# print '====='
# ave_ref = np.unique(gains_cal['a2'])
# for a in ave_ref:
#     print a, np.sum(gains_cal['a2']==a)/num_antennas


make_plot = True
make_plot_1 = False
if make_plot:
    if make_plot_1:
        fig, axes = pp.subplots(6, 1, sharex=True, sharey=False, figsize=(12,10))

        # -------------------
        ax = axes[0]
        ax.set_title('Gains for antenna %i' % a1)
        x = gains_cor_a['time']-gains_cor_a['time'][0]
        y = np.angle(gains_cor_a['cparam'][:, 0])*(180./np.pi)
        y = -y
        ax.plot(x, y, 'ro-', markerfacecolor='None', markeredgecolor='r', ms=3, lw=2,
                label='corruption')
        x = gains_cal_a['time']-gains_cal_a['time'][0]
        y = np.angle(gains_cal_a['cparam'][:, 0])*(180./np.pi)
        ax.plot(x, y, 'b+--', markerfacecolor='None', markeredgecolor='b', ms=3, lw=2,
                label='calibrated')
        # x = gains_cal_ave_a['time']-gains_cal_ave_a['time'][0]
        # y = np.angle(gains_cal_ave_a['cparam'][:, 0])*(180./np.pi)
        # ax.plot(x, y, 'g.', markerfacecolor='None', markeredgecolor='g', ms=3, lw=2,
        #         label='BDA calibrated')
        ax.set_ylabel('gain phase [deg]')
        ax.legend(ncol=2, prop={'size':8})
        ax.grid()


        # -------------------
        ax = axes[1]
        x = gains_cor_a['time']-gains_cor_a['time'][0]
        y = 1./np.abs(gains_cor_a['cparam'][:, 0])
        ax.plot(x, y, 'ro-', markerfacecolor='None', markeredgecolor='r', ms=3, lw=2,
                label='corruption')
        x = gains_cal_a['time']-gains_cal_a['time'][0]
        y = np.abs(gains_cal_a['cparam'][:, 0])
        ax.plot(x, y, 'b+--', markerfacecolor='None', markeredgecolor='b', ms=3, lw=2,
                label='calibrated')
        # x = gains_cal_ave_a['time']-gains_cal_ave_a['time'][0]
        # y = np.abs(gains_cal_ave_a['cparam'][:, 0])
        # ax.plot(x, y, 'g.', markerfacecolor='None', markeredgecolor='g', ms=3, lw=2,
        #         label='BDA calibrated')
        ax.set_ylabel('gain amplitude')
        ax.legend(ncol=2, prop={'size':8})
        ax.grid()

        # -------------------
        ax = axes[2]
        x = gains_cal_a['time']-gains_cal_a['time'][0]
        y = np.abs(gains_cal_a['snr'][:, 0])
        ax.plot(x, y, 'bo', markerfacecolor='None', markeredgecolor='b', ms=3, lw=1,
                label='calibrated')
        # x = gains_cal_ave_a['time']-gains_cal_ave_a['time'][0]
        # y = np.abs(gains_cal_ave_a['snr'][:, 0])
        # ax.plot(x, y, 'g.', markerfacecolor='None', markeredgecolor='g', ms=3, lw=2,
        #         label='BDA calibrated')
        ax.set_ylabel('gain snr')
        ax.legend(ncol=2, prop={'size':8})
        ax.grid()

        # -------------------
        ax = axes[3]
        x = gains_cal_a['time']-gains_cal_a['time'][0]
        y = np.abs(gains_cal_a['paramerr'][:, 0])
        ax.plot(x, y, 'bo', markerfacecolor='None', markeredgecolor='b', ms=3, lw=2,
                label='calibrated')
        # x = gains_cal_ave_a['time']-gains_cal_ave_a['time'][0]
        # y = np.abs(gains_cal_ave_a['paramerr'][:, 0])
        # ax.plot(x, y, 'g.', markerfacecolor='None', markeredgecolor='g', ms=3, lw=2,
        #         label='BDA calibrated')
        ax.set_ylabel('gain paramerr')
        ax.legend(ncol=2, prop={'size':8})
        ax.grid()

        # -------------------
        ax = axes[4]
        x = gains_cal_a['time']-gains_cal_a['time'][0]
        y = np.abs(gain_check)
        ax.plot(x, y, 'bo', markerfacecolor='None', markeredgecolor='b', ms=3, lw=2,
                label='check')
        # x = gains_cal_ave_a['time']-gains_cal_ave_a['time'][0]
        # y = np.abs(gains_cal_ave_a['paramerr'][:, 0])
        # ax.plot(x, y, 'g.', markerfacecolor='None', markeredgecolor='g', ms=3, lw=2,
        #         label='BDA calibrated')
        ax.set_ylabel('gain amp check')
        ax.grid()

        # -------------------
        ax = axes[5]
        x = gains_cal_a['time']-gains_cal_a['time'][0]
        y = np.angle(gain_check)*(180./np.pi)
        ax.plot(x, y, 'bo', markerfacecolor='None', markeredgecolor='b', ms=3, lw=2,
                label='check')
        # x = gains_cal_ave_a['time']-gains_cal_ave_a['time'][0]
        # y = np.abs(gains_cal_ave_a['paramerr'][:, 0])
        # ax.plot(x, y, 'g.', markerfacecolor='None', markeredgecolor='g', ms=3, lw=2,
        #         label='BDA calibrated')
        ax.set_ylabel('gain phase check [deg]')
        ax.grid()


    fig = pp.figure()



    for i in range(0, 1):
        # Obtain gains for antenna a1 with phase referene of antenna a2
        a1 = i
        a2 = 0
        gains_cal_ave_a = gains_cal_ave[gains_cal_ave['a1']==a1]
        gains_cal_ave_a = gains_cal_ave_a[gains_cal_ave_a['a2']==a2]
        a2 = 0
        gains_cal_a = gains_cal[gains_cal['a1']==a1]
        gains_cal_a = gains_cal_a[gains_cal_a['a2']==a2]
        a2 = -1
        gains_cor_a = gains_cor[gains_cor['a1']==a1]
        gains_cor_a = gains_cor_a[gains_cor_a['a2']==a2]

        # Remove flagged antenna
        gains_cal_ave_a = gains_cal_ave_a[gains_cal_ave_a['flag'][:,0]==0]
        gains_cal_a = gains_cal_a[gains_cal_a['flag'][:,0]==0]

        pp.clf()
        ax = fig.add_subplot(111)
        ax.set_title('Gains antenna %i' % a1)
        x = gains_cal_ave_a['time']-gains_cal_ave_a['time'][0]
        y = np.abs(gains_cal_ave_a['cparam'][:, 0])
        ax.plot(x, y, 'g.', markerfacecolor='None', markeredgecolor='g', ms=3, lw=2,
                label='BDA calibrated')
        x = gains_cal_a['time']-gains_cal_a['time'][0]
        y = np.abs(gains_cal_a['cparam'][:, 0])
        ax.plot(x, y, 'b.', markerfacecolor='None', markeredgecolor='b', ms=3, lw=2,
                label='calibrated')
        ax.legend(ncol=2, prop={'size':8})
        ax.grid()
        pp.draw()
        pp.show()
        time.sleep(0.2)
