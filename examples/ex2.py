#!/usr/bin/env python
#
# Purpose:
#     Solving Helmholtz eq. on a square domain using Legendre polynomial collocation.
# Problem being solved: 
#    Helmholtz eq. u_xx + u_yy + (k^2)u = f, on [-1,1]x[-1,1]
# Boundary conditions: 
#     Homogenous Dirichlet BC's
#
from pysplib import *
import numpy as np
from scipy.linalg import solve
from scipy.interpolate import interp2d
from matplotlib import pyplot as plt

import pylab as p
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm

# polynomial order
n=32

# Legendre Gauss-Lobatto nodes
x,vn=zelegl(n)

# Derivative matrix for Legendre at Gauss-Lobatto nodes
D=dmlegl(n,x,vn,n)

# Tensor product grid
y=x 
xx,yy = np.meshgrid(x[1:n], y[1:n])
xx=xx.flatten(1)
yy=yy.flatten(1)

# Construct rhs vector
f=np.exp(-10.*((yy-1.)**2+(xx-0.5)**2))

# Construct differential operator
D2=np.dot(D,D)
D2=D2[1:n,1:n]
I=np.eye(n-1)
k=9.
L=np.kron(I,D2)+np.kron(D2,I)+k**2*np.eye((n-1)**2)# Helmholtz operator

# Solve the system
u=solve(L,f)

# Reshape long 1D results to 2D grid:
uu=np.zeros((n+1,n+1))
uu[1:n,1:n] = u.reshape(n-1,n-1).transpose()
xx,yy = np.meshgrid(x,y)
value = uu[n/2,n/2]

# Interpolate to finer grid and plot
uu = uu.flatten(1)
uuu = interp2d(x,y,uu,kind='cubic')

a = np.linspace(-1.,1.,100)
b = np.linspace(-1.,1.,100)
xxx,yyy = np.meshgrid(a,b)
newfunc = uuu(a,b)

# Plot results
# wireframe:
fig=p.figure()
ax = p3.Axes3D(fig)
ax.plot_wireframe(xxx,yyy,newfunc)
ax.text(-0.8,0.5,.022,'u(0,0) = %13.11f' % value)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
p.show()

# mlab plotting
from enthought.mayavi import mlab
mlab.surf(a,b,newfunc)
mlab.show()
