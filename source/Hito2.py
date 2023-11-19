
from numpy import array, zeros, linspace, arange
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton
from Schemes.Schemes2 import *

## HITO 2

def F_Cauchy(U, t=0): # df/dx = x

    x = U[0]
    return x
def ODE(U, t=0): # df/dx = x
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    vx=1-4*x+x^2*y^2
    vy=3*x-x^2*y

    returnarray( [vx, vy] )

def F_Kepler(U, t=0):
 
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )


def Euler(F,U,t0, tf):
    dt = tf-t0
    return U + dt * F(U,t0)

def RK4(F,U,t0, tf ):
    dt = tf-t0
    k1 = F( U, t0 )
    k2 = F( U +k1*dt/2, t0 + dt/2 )
    k3 = F( U +k2*dt/2, t0 + dt/2 )
    k4 = F( U + k3*dt, t0 + dt )
     
    k = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    return U + dt * k

def Crank_Nicolson(F, U, t0, tf):
    dt = tf-t0
    def g(x):
        return x - a -dt/2 * (F(a,t0) + F(x,tf))
    a = U
    
    return newton(g, U)

def Inverse_Euler(F, U, t0, tf):
    dt = tf-t0
    def g(x):
        return x - a - dt * F(x,tf) 
    a = U
    
    return newton(g, U)

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