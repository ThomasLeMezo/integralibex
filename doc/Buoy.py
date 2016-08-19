#!/usr/bin/python
#
"""
Created on Sun Nov 15 16:22:16 2015

@author: lemezoth
"""

import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

def vanDerPol(y, t):
    dy0 = y[1]
    dy1 = 1.0*(1-y[0]**2)*y[1]-y[0]
    return [dy0, dy1]

def chemin(x, y, color):
    t = np.linspace(0, 100, 10000)
    r = odeint(vanDerPol, (x, y), t)
    plt.plot(r[:,0], r[:,1], color)

def integration(X, dt, tmax, color):
	nb_item = int(round(tmax/dt))
	r = np.zeros([2, nb_item])
	for i in range(0, nb_item):
		[dx0, dx1] = vanDerPol(X, 0)
		X[0] += dt*dx0
		X[1] += dt*dx1
		r[:,i] = np.array([X[0], X[1]])
	plt.plot(r[0,:], r[1,:], color)


# chemin(0.00001, 0.00001, 'k')
# chemin(0.00001, -0.00001, 'k')
# chemin(-0.000001, 0.00001, 'k')
# chemin(-0.00001, -0.00001, 'k')

# chemin(5.0, 5.0, 'k:')
# chemin(0, -4, 'k:')
integration([3.0, -3.0], 0.001, 10, 'k')

plt.axis('equal')
plt.xlim((-4, 4))
plt.ylim((-4, 4))
plt.show()