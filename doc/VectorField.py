#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from numpy import *

# Set limits and number of points in grid

##########  Part 1 ##########
phi, d = np.meshgrid(np.arange(-np.pi, np.pi, .2), np.arange(0.1, 10, .2))
U = np.where(np.cos(phi) <= 1/np.sqrt(2.0), np.sin(phi)/d+1.0, (1.0/d-1.0)*np.sin(phi))
V = -np.cos(phi)
plt.figure()
# Q = plt.quiver(U, V)

# plt.show()

########## Part 2 ##########
# t = np.linspace(0.0, 2*np.pi, 1000)
# # for c in drange(0.0, 1.0, 0.1):
# c=0.78
# d = np.sqrt((cos(t)+c)**2+(sin(t))**2)
# phi = arctan2(cos(t), -sin(t))+np.pi-arctan2(sin(t), cos(t)+c)
# phi_mod = 2*arctan(tan(phi/2.0))

# plt.plot(phi_mod, d)

# plt.title('c=')

# l, r, b, t = plt.axis()
# dx, dy = r - l, t - b
# plt.axis([l - 0.05*dx, r + 0.05*dx, b - 0.05*dy, t + 0.05*dy])
# plt.show()

########## Part 3 ##########
c=0.705
t = np.linspace(0.0, 2*np.pi, 1000)
phi = arctan2(cos(t), -sin(t))+np.pi-arctan2(sin(t), cos(t)+c)
phi_mod = 2*arctan(tan(phi/2.0))
plt.plot(t, phi_mod)

limit = -np.ones(size(t))*pi/4.0
plt.plot(t, limit, 'r')

plt.show()