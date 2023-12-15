
from numpy import array, zeros, linspace, arange
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton
from Temp_Schemes.Schemes import *
from ODEs.ODE import *

## HITO 2

def Integrate_ODE(U, F, t, scheme):

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

N = 1000
dt = 0.01                        
t = linspace(0, N*dt, N+1)

#Keplerian Orbit initial conditions

#Initial conditions Cauchy problem
U0_c = array(zeros((2,len(t)-1)))
U0_c[:,0] = array( [1, 0] )

#Initial conditions Keplerian orbit problem, circular orbit
U0_k = array(zeros((4,len(t)-1)))
U0_k[:,0] = array( [1, 0, 0, 1] )

U_c = Integrate_ODE(U0_c, ODE, t, RK4)
U_k = Integrate_ODE(U0_k, F_Kepler, t, Inverse_Euler)

#Plot
plt.plot(U_c[0,:],U_c[1,:])
plt.axis('equal')
plt.show()
plt.plot(U_k[0,:],U_k[1,:])
plt.axis('equal')
plt.show()