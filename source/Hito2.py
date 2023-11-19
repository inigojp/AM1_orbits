
from numpy import array, zeros, linspace, arange
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton
from Temp_Schemes.Schemes import *
from ODEs.ODE import F_Kepler

## HITO 2

def Integrate_ODE(U, F, t, scheme):

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

N = 1000
dt = 0.01                        
t = linspace(0, N*dt, N+1)

#Keplerian Orbit initial conditions

U = array(zeros((4,len(t)-1)))
U[:,0] = array( [1, 0, 0, 1] )

U = Integrate_ODE(U, F_Kepler, t, RK4)


#Plot
plt.plot(U[0,:],U[1,:])
plt.axis('equal')
plt.show()