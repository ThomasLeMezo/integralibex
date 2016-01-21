#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma

# Set limits and number of points in grid

phi, d = np.meshgrid(np.arange(-np.pi, np.pi, .2), np.arange(0.1, 10, .2))
U = np.where(np.cos(phi) >= np.sqrt(2.0)/2.0, np.sin(phi)/d+1.0, (1.0/d-1.0)*np.sin(phi))
V = -np.cos(phi)

# 1
plt.figure()
Q = plt.quiver(U, V)

l, r, b, t = plt.axis()
dx, dy = r - l, t - b
plt.axis([l - 0.05*dx, r + 0.05*dx, b - 0.05*dy, t + 0.05*dy])

plt.title('Vector Field')
plt.show()