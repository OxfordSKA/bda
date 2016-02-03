# coding=utf-8

import numpy as num
import numpy.linalg as LA


def stefcal1(a, b, tol, niter, gstart=None):
    """
    Python implementation of StefCal calibration alogorithm

    Description:
        Minimize  || a - G' * b * G ||_F where G = diag(g)

    Args:
        a (array_like): Observed visibilities (full matrix)
        b (array_like): Model visibilities (full matrix)
        tol (float): Required tolerance (eg. 1.0e-8)
        niter (int): Maximum number of iterations (eg. 50)
        gstart (array_like): Initial values for g. If not present, set g = 1
                             Useful as a starting point if, for example, gains
                             at a previous time are known.

    Returns:
        g[], nit, dg
        g (array_like): computed gains, g(i, 1) corresponds to the ith receiver.
        nit (int): number of iterations required.
        dg (float): convergence achieved.
    """

    n = a.shape[0]

    if gstart is not None:
        g = gstart
    else:
        g = num.ones(n, num.complex128)

    dg = 1.0e30
    nit = niter
    omega = 0.5
    f0 = 1 - omega
    f1 = omega

    for iter in range(niter):
        gold = num.copy(g)
        nit = niter
        if iter < 2:
            for j in range(n):
                z = num.conj(gold)*b[:, j]
                w = num.dot(num.conj(z), z)
                t = num.dot(num.conj(z), a[:, j])
                t = t/w
                g[j] = t
            dg = LA.norm(g-gold) / LA.norm(g)
            if dg <= tol:
                nit = iter
                break
        else:
            for j in range(n):
                z = num.conj(gold)*b[:, j]
                w = num.dot(num.conj(z), z)
                t = num.dot(num.conj(z), a[:, j])
                t = t/w
                g[j] = t
            if iter % 2 == 1 and iter > 0:
                dg = LA.norm(g - gold) / LA.norm(g)
                if dg <= tol:
                    nit = iter
                    break
                else:
                    g = f0 * g + f1 * gold

    p = num.conj(g[0]) / num.abs(g[0])
    g = p * g

    return g, nit, dg
