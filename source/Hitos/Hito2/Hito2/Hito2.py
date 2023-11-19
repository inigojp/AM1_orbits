
from numpy import array, zeros, linspace, arange
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton

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

def Euler(F,i):
    
    return U[:,i] + dt * F(U[:,i])
    
def Crank_Nicolson(F,i):
    
    def g(x):
        return x - a -dt/2 * (F(a) + F(x))
    a = U[:,i]
    
    return newton(g, U[:,i])

def RK4(F,i):
    
    k1 = F( U[:,i], t )
    k2 = F( U[:,i] +k1*dt/2, t + dt/2 )
    k3 = F( U[:,i] +k2*dt/2, t + dt/2 )
    k4 = F( U[:,i] + k3*dt, t + dt )
     
    k = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    return U[:,i] + dt * k

def Inverse_Euler(F,i):
    def g(x):
        return x - a - dt * F(x,t) 
    a = U[:,i]
    
    return newton(g, U[:,i])

def Integrate_ODE(U, F, scheme):
    
    for i in range(0, N-1):
        U[:,i+1] = scheme(F,i)
        
    return U

N = 100
dt = 0.01                        
t = arange(N)

#Keplerian Orbit initial conditions
# U = array(zeros((4,N)))
# U[:,0] = array( [1, 0, 0, 1] )

#Cauchy df/dx = x initial conditions
U = array(zeros((1,N)))
U[:,0] = 0.2

U = Integrate_ODE(U, F_Cauchy, RK4)

plt.plot(U[0,:])

plt.show()