# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# x = np.arange(0, 2*np.pi + np.pi/4, 2 * np.pi / 8)
# y = np.sin(x)
# tck = interpolate.splrep(x, y, s=0)
# xnew = np.arange(0, 2*np.pi, np.pi / 50)
# ynew = interpolate.splev(xnew, tck, der=0)
#
#
# plt.figure()
# plt.plot(x, y, 'bx-')
# plt.plot(xnew, ynew, 'g', lw=3.0)
# plt.plot(xnew, np.sin(xnew), 'r--')
# plt.plot(x, y, 'b')
# plt.legend(['Linear', 'Cubic Spline', 'True'])
# plt.axis([-0.05, 6.33, -1.05, 1.05])
# plt.title('Cubic-spline interpolation')
# plt.show()

x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
y = np.sin(x)
s = interpolate.InterpolatedUnivariateSpline(x, y)
xnew = np.arange(0, 2*np.pi, np.pi/50)
ynew = s(xnew)
print s.integral(0, np.pi)
print s.integral(np.pi, 2*np.pi)
print s.integral(0, 2.0*np.pi)

plt.figure()
plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
plt.legend(['Linear', 'InterpolatedUnivariateSpline', 'True'])
plt.axis([-0.05, 6.33, -1.05, 1.05])
plt.title('InterpolatedUnivariateSpline')
plt.show()
