# -*- coding: utf-8 -*-
"""Module to plot a gains pickle file."""


import numpy as np
import matplotlib.pyplot as plt
import os
import pickle


if __name__ == "__main__":

    station = 80

    gains_pickle1 = os.path.join('out_default', 'vis', 'gains.pickle')
    gains1 = pickle.load(open(gains_pickle1))
    plt.plot(np.real(gains1[station]), 'r+-')
    plt.plot(np.imag(gains1[station]), 'r+-')

    gains_pickle2 = os.path.join('out_new', 'vis', 'gains.pickle')
    gains2 = pickle.load(open(gains_pickle2))
    plt.plot(np.real(gains2[station]), 'bx--')
    plt.plot(np.imag(gains2[station]), 'bx--')

    plt.show()

    print gains1[station][0]
    print gains2[station][0]






