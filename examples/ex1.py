#!/usr/bin/env python
#
# Purpose:
#     Solving two point BVP using Legendre polynomial collocation.
# Problem being solved: 
#    -d2U/dX2 + tau*U = f,  for x in [-1,1]
# Dirichlet boundary conditions: 
#     U(-1) = sig1; U(1) = sig2
# Here we set f=1., tau=1., sig1=sig2=0.
#
from pysplib import *
import numpy as np
from scipy.linalg import solve
from matplotlib import pyplot as plt

# Constants
sig1=0.
sig2=0.
tau=1.
ef=1.
exp=np.exp(1)
cost=exp/(1.+exp*exp)

# polynomial order
n=10

# Legendre Gauss-Lobatto nodes
et,vn=zelegl(n)

# ...and weights
wt=welegl(et,vn,n)

# Derivative matrix for Legendre at Gauss-Lobatto nodes
dma=dmlegl(n,et,vn,n)

# Construct differential operator
l=-np.dot(dma,dma)+tau*np.eye(n+1)

# Change the entries according to BCs
l[0,:]=0.
l[n,:]=0.
l[0,0]=1.
l[n,n]=1.

# Construction of the rhs vector
fu=np.zeros(n+1)
fu[:]=ef
fu[0]=sig1
fu[n]=sig2

# Solve the system
u=np.zeros(n+1)
u=solve(l,fu)

# Errr vector
ei=np.exp(et)
qn=u-(1.-cost*(ei+1./ei))

# Norm of the error vector
qi,qs,qm = nolegl(vn,qn,wt,n=(len(vn)-1))

# Interpolation to finer grid for plotting
xx=np.linspace(-1.,1.,100)
uu = np.polyval(np.polyfit(et,u,n),xx)

# Plot the result
plt.title('integral err = %e' % qi)
plt.plot(xx,uu,'b')
plt.show()

