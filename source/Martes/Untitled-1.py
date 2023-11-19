
from numpy import array, zeros, linspace, arange
from Temporal_schemes.Schemes2 import Inverse_Euler
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton
import math
## Integración sistemas dinámicos

## ODEs
# Sistema caótico con k<1 
def ODE_DUffing(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=-x**3-0.01*y+math.cos(t)

    return array( [vx, vy] )

def ODE_Van_Der_Por_1(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=-0.2*(x**2-1)*y-x

    return array( [vx, vy] )

def ODE_Van_Der_Por_2(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=-5*(x**2-1)*y-x

    return array( [vx, vy] )

def ODE_Rayleigh(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=x-x**3-y+0.6*math.cos(t)

    return array( [vx, vy] )

def ODE_Rayleigh2(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=x-x**3-y+3*math.cos(t)

    return array( [vx, vy] )

## Esquemas numéricos


##Integrar 
def Integrate_ODE(U, F, t, scheme):

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

N = 30000
d_t = 0.01                        
t = linspace(0,N*d_t, N+1)

#Keplerian Orbit initial conditions
U = array(zeros((2,len(t)-1)))
U[:,0] = array( [1, 0.01] )
#U_0 = array( [1.1, 3] )
U2 = array(zeros((2,len(t)-1)))
U2[:,0] = array( [1, 0.01] )
#Cauchy df/dx = x initial conditions
#U = array(zeros((1,N)))
#U[:,0] = 1.1

#U = Integrate_ODE(U, ODE_Rayleigh, t, RK4)
U2 = Integrate_ODE(U2, ODE_Rayleigh2, t, Inverse_Euler)


f1 = plt.figure()
#plt.plot(U[0,5000:N],U[1,5000:N])
plt.plot(U2[0,10000:N],U2[1,10000:N])
#f2 = plt.figure()
#plt.plot(U2[0,0000:N])
#plt.plot(U2[1,5000:N])

#plt.plot(U2[0,:])

plt.show()