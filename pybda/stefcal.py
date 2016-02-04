# coding=utf-8
"""
Python implementation of StefCal calibration algorithm.
"""


import numpy
import numpy.linalg


def stefcal1(a, b, tol=1.0e-8, niter=50, gstart=None):
    """
    Python implementation of StefCal calibration algorithm.

    Description:
        Minimize
            ``|| a - G' * b * G ||_F``
        where
            ``G = diag(g)``

    Args:
        a (array_like): Observed visibilities (full matrix)
        b (array_like): Model visibilities (full matrix)
        tol (float): Required tolerance (eg. 1.0e-8)
        niter (int): Maximum number of iterations (eg. 50)
        gstart (array_like): Initial values for g. If not present, set g = 1
                             Useful as a starting point if, for example, gains
                             at a previous time are known.

    Returns:
        tuple of g, nit, dg

        g (array_like): computed gains, g[i] corresponds to the ith receiver.
        nit (int): number of iterations required.
        dg (float): convergence achieved.
    """
    a = numpy.asarray(a)
    b = numpy.asarray(b)
    n = a.shape[0]

    if gstart is not None:
        g = numpy.asarray(gstart)
    else:
        g = numpy.ones(n, numpy.complex128)

    dg = 1.0e30
    nit = niter
    omega = 0.5
    f0 = 1 - omega
    f1 = omega

    for i in range(niter):
        g_old = numpy.copy(g)
        for j in range(n):
            z = numpy.conj(g_old) * b[:, j]
            g[j] = numpy.dot(numpy.conj(z), a[:, j]) / \
                numpy.dot(numpy.conj(z), z)
        if i < 2:
            dg = numpy.linalg.norm(g - g_old) / numpy.linalg.norm(g)
            if dg <= tol:
                nit = i
                break
        else if i % 2 == 1:
            dg = numpy.linalg.norm(g - g_old) / numpy.linalg.norm(g)
            if dg <= tol:
                nit = i
                break
            else:
                g = f0 * g + f1 * g_old

    p = numpy.conj(g[0]) / numpy.abs(g[0])
    g = p * g

    return g, nit, dg
