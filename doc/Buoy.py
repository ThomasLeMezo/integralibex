# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 16:22:16 2015

@author: lemezoth
"""

import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

def vanDerPol(y, t):
    mu = 1.0
    dy0 = y[1]
    dy1 = mu*(1-y[0]**2)*y[1]-y[0]
    return [dy0, dy1]

def chemin(x, y, color):
    t = np.linspace(0, 100, 10000)
    r = odeint(vanDerPol, (x, y), t)
    plt.plot(r[:,0], r[:,1], color)


chemin(0.00001, 0.00001, 'k')
chemin(0.00001, -0.00001, 'k')
chemin(-0.000001, 0.00001, 'k')
chemin(-0.00001, -0.00001, 'k')

chemin(5.0, 5.0, 'k:')
chemin(0, -4, 'k:')

plt.axis('equal')
plt.xlim((-4, 4))
plt.ylim((-4, 4))
plt.show()