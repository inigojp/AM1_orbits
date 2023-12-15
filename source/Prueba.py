from numpy import array, zeros, log10, ones, vstack  
from numpy.linalg import norm, lstsq
from Temp_Schemes.Schemes import *
from ODEs.ODE import F_Kepler
import matplotlib.pyplot as plt

def Integrate_ODE(U0, F, t, scheme):
    N, Nv=  len(t)-1, len(U0)
    U = zeros( (Nv, N), dtype=float)
    U[:,0] = U0

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

t0 = 0
tf = 1000
dt = 0.01         
Np = int((tf-t0)/dt)               
t = linspace(t0, tf, Np)


U0_k = array( [1, 0, 0, 1] )

U_k = Integrate_ODE(U0_k, F_Kepler, t, RK4)

plt.plot(U_k[0,:],U_k[1,:])
plt.axis('equal')
plt.show()