# -*- coding: utf-8 -*-

import os
import numpy as np
import pickle
import progressbar
import matplotlib.pyplot as plt

def load_ms(ms):
    dtype = [('a1','i4'), ('a2', 'i4'), ('data', 'c16'), ('t', 'i4')]
    tb.open(ms)
    num_rows = tb.nrows()
    values = np.zeros((num_rows,), dtype=dtype)
    values['data'] = np.squeeze(tb.getcol('DATA'))
    values['a1'] = np.squeeze(tb.getcol('ANTENNA1'))
    values['a2'] = np.squeeze(tb.getcol('ANTENNA2'))
    time = tb.getcol('TIME')
    time -= time[0]
    time /= 0.08
    time_index = np.array(np.round(time), dtype='i4')
    values['t'] = time_index
    tb.close()
    return values


model_ms = os.path.join('out_default', 'vis', 'model.ms')
default_ms = os.path.join('out_default', 'vis', 'corrupted.ms')
new_ms = os.path.join('out_new', 'vis', 'corrupted.ms')
gains = os.path.join('out_default', 'vis', 'gains.pickle')

model = load_ms(model_ms)
default = load_ms(default_ms)
new = load_ms(new_ms)
g = pickle.load(open(gains))

dtype = [('diff', 'c16'), ('abs_diff', 'f8')]
dtype = dtype + default.dtype.descr
vis = np.empty(default.shape, dtype=dtype)
for name in default.dtype.names:
    vis[name] = default[name]

diff = default['data'] - new['data']
vis['diff'] = diff
vis['abs_diff'] = np.abs(diff)

print '*' * 80
print 'abs diff :', vis['abs_diff'][0]
print 'diff     :', vis['diff'][0]
print 'data     :', vis['data'][0]
print 'a1       :', vis['a1'][0]
print 'a2       :', vis['a2'][0]
print 't        :', vis['t'][0]
vis.sort(order='abs_diff')
vis = vis[::-1]
print ''
print 'abs diff :', vis['abs_diff'][0]
print 'diff     :', vis['diff'][0]
print 'data     :', vis['data'][0]
print 'a1       :', vis['a1'][0]
print 'a2       :', vis['a2'][0]
print 't        :', vis['t'][0]
print '*' * 80


p = 86
q = 101
t = 141
model_pq = model[new['a1'] == p]
model_pq = model_pq[model_pq['a2'] == q]
new_pq = new[new['a1'] == p]
new_pq = new_pq[new_pq['a2'] == q]
default_pq = default[default['a1'] == p]
default_pq = default_pq[default_pq['a2'] == q]
# vis_pq = vis[vis['a1'] == p]
# vis_pq = vis_pq[vis_pq['a2'] == q]

print 'model   ', model_pq['data'][t]
print 'default ', default_pq['data'][t]
print 'new     ', new_pq['data'][t]

gp = g[p][t]
gq = g[q][t]
print '-' * 80
print '%3i' % p, gp
print '%3i' % q, gq
Vpq = 1.0/gp * model_pq['data'][t] * np.conj(1.0/gq)
print 'Vpq', Vpq
print '-' * 80


# plt.plot(g[p],'+--')
# plt.plot(g[q],'x--')
# plt.show()

#
#
# # print '-' * 80
# # mpq = model['data'][0]
# # gp = g[model['a1'][0]][0]
# # gq = g[model['a2'][0]][0]
# # # gp = np.complex64(gp)
# # # gq = np.complex64(gq)
# # vpq = 1./gp * mpq * np.conj(1./gq)
# # print mpq, gp, gq
# # print vpq
# # print 'diff1:', vpq - default['data'][0]
# # print 'diff2:', vpq - new['data'][0]
# # print '-' * 80
# #
# #
# #
# # plt.plot(vis_pq['abs_diff'], '+')
# # # plt.plot(np.abs(new_pq['data']), '+')
# # # plt.plot(np.abs(default_pq['data']), 'x')
# # plt.show()
# #
# # diff_data = default['data'] - new['data']
# # diff_a1 = default['a1'] - new['a1']
# # diff_a2 = default['a2'] - new['a2']
# # print np.max(diff_a1)
# # print np.max(diff_a2)
# # print np.max(np.abs(diff_data))
