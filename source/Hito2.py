
from numpy import array, zeros, linspace, arange
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton
from Temp_Schemes.Schemes import *
from ODEs.ODE import *

## HITO 2

def Integrate_ODE(U, F, t, scheme):

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i], t[i+1] - t[i], t[i])
        
    return U

#Inputs
t0 = 0
tf = 100
dt = 0.01

# Discretización de tiempo
N = int((tf-t0) / dt)                    
t = linspace(t0, tf, N)  # Tiempo inicial, tiempo final, número de puntos

#Keplerian Orbit initial conditions

#Initial conditions Cauchy problem
U0_c = array(zeros((2,len(t))))
U0_c[:,0] = array( [1, 0] )
U0_r = array(zeros((2,len(t))))
U0_r[:,0] = array( [1, 0] )



U_c = Integrate_ODE(U0_c, ODE, t, RK4)
U_r = Integrate_ODE(U0_r, ODE_Rayleigh, t, Inverse_Euler)

#Plot
plt.plot(U_c[0,:],U_c[1,:])
plt.axis('equal')
plt.show()
plt.plot(U_r[0,:],U_r[1,:])
plt.axis('equal')
plt.show()