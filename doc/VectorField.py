#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from numpy import *

# Set limits and number of points in grid

##########  Part 1 ##########
# phi, d = np.meshgrid(np.arange(-np.pi, np.pi, .2), np.arange(0.1, 10, .2))
# U = np.where(np.cos(phi) <= 1/np.sqrt(2.0), np.sin(phi)/d+1.0, (1.0/d-1.0)*np.sin(phi))
# V = -np.cos(phi)
# plt.figure()
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
# c=1.0/sqrt(2.0)
# t = np.linspace(0.0, 2*np.pi, 1000)
# phi = arctan2(cos(t), -sin(t))+np.pi-arctan2(sin(t), cos(t)+c)
# phi_mod = 2*arctan(tan(phi/2.0))
# plt.plot(t, phi_mod)

# limit = -np.ones(size(t))*pi/4.0
# plt.plot(t, limit, 'r')

# plt.show()

#### TEST ####

#### CIRCLE
# x1, x2 = np.meshgrid(np.arange(-2.0,2.0, .1), np.arange(-2.0,2.0, .1))
# U = x1-(x1+x2)*(x1**2+x2**2)
# V = x2+(x1-x2)*(x1**2+x2**2)
# plt.figure()
# Q = plt.quiver(U, V)
# plt.show()

#### Parrilo
# x1, x2 = np.meshgrid(np.arange(-200.0,200.0, 10.0), np.arange(-200.0,200.0,10.0))
# U = -x1+(1+x1)*x2
# V = -(1+x1)*x1

#### Ratschan 3
# x1, x2 = np.meshgrid(np.arange(-200.0,200.0, 10.0), np.arange(-200.0,200.0,10.0))
# x1, x2 = np.meshgrid(np.arange(-1.0,2.0, 0.1), np.arange(-1.0,1.0,0.1))
# U = -4*x1*x1*x1+6*x1*x1-2*x1
# V = -2*x2

#### Genesio
#x1, x2 = np.meshgrid(np.arange(-70,70,1.0), np.arange(-1000,1000,100))
#x1, x2 = np.meshgrid(np.arange(-20,10,0.5), np.arange(-15,5,0.5))
# x1, x2 = np.meshgrid(np.arange(-0.5,0.5,0.05), np.arange(-0.5,0.5,0.05))
#U=(-x1+x2)
#V=(0.1*x1-2*x2-x1*x1-0.1*x1*x1*x1)

#### Magnet ?
# x1, x2 = np.meshgrid(np.arange(-5.0,5.0,0.1), np.arange(-5.0,5.0,0.1))
# U = x1
# V = 1.0/(x1-1.0)-1.0/(x1+1.0)

##### Waves
x1, x2 = np.meshgrid(np.arange(-10,10,0.5), np.arange(0.0,10.0,0.5))
U = x1+exp(x2)*sin(x1+2)
V = x2-exp(x2)*cos(x1+2)


#### Van Der Pol
# x1, x2 = np.meshgrid(np.arange(-1.0,13.0, .3), np.arange(-6.0,6.0, .2))
# x1, x2 = np.meshgrid(np.arange(-4.0,4.0, .4), np.arange(-4.0,4.0, .4))
# U = x2
# V = (1-x1**2)*x2-x1

#### Car on the hill ?
# U = x2
# V = -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2

#### Host & parasites
# x1, x2 = np.meshgrid(np.arange(-1.0,2.0, 0.1), np.arange(-1.0,2.0, 0.1))
# x1, x2 = np.meshgrid(np.arange(-1.0,2.0, 0.1), np.arange(-1.0,2.0, 0.1))
# U = (1-x2)*x1
# V = x2*(1-2*x2/(1+x1))

#### Drawing

coeff = np.maximum(np.minimum(V**2+U**2, 1.5*np.ones(V.shape)), 0.5*np.ones(V.shape))

U = U / (sqrt(U*U+V*V))
V = V / (sqrt(U*U+V*V))

# U = coeff * U / (sqrt(U*U+V*V))
# V = coeff * V / (sqrt(U*U+V*V))


# plt.figure()
Q = plt.quiver(x1, x2, U, V, width=1e-3)
plt.show()
