from numpy import array, zeros, linspace
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton

N = 1000                         
dt = 0.01                        
t = array( zeros(N) )

def F_Kepler(U,t):
 
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )

#EULER
U = array( [1, 0, 0, 1] )
x = array( zeros(N) )
y = array( zeros(N) )
x[0] = U[0]
y[0] = U[1]

for i in range(0, N):
    F = F_Kepler( U,t )
    U = U + dt * F
    x[i] = U[0]
    y[i] = U[1]

plt.axis('equal')
plt.plot(x, y)
plt.show()

# CRANK-NICOLSON
U = array(zeros((4,N)))
U[:,0] = array( [1, 0, 0, 1] )

def g(x):
    return x - a -dt/2 * (F_Kepler(a,t) + F_Kepler(x,t))

for n in range (0, N-1):
    a = U[:,n]
    U[:,n+1] = newton(g, U[:,n])
    
plt.axis('equal')
plt.plot(U[0,:],U[1,:])
plt.show()

# #RUNGE-KUTTA 4TH ORDER
U = array( [1, 0, 0, 1] )
x_rk = array(zeros(N))
y_rk = array(zeros(N))
x_rk[0] = U[0]
y_rk[0] = U[1]

for i in range(0,N):
    

    k1 = F_Kepler( U, t)
    k2 = F_Kepler( U +k1*dt/2, t + dt/2 )
    k3 = F_Kepler( U +k2*dt/2, t + dt/2)
    k4 = F_Kepler( U + k3*dt, t + dt )
     
    k = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    U = U + dt * k
    
    x_rk[i] = U[0]
    y_rk[i] = U[1]
    
plt.axis('equal')
plt.plot(x_rk, y_rk)
plt.show()