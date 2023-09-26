from numpy import array, zeros, linspace
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

N = 10000                       # PASOS DE INTEGRACION
dt = 0.01                        # INTERVALO ENTRE PASOS
t_values = array(zeros(N))

def F_Kepler(U):
 
    x, y, vx, vy = [U[0,:], U[1,:], U[2,:], U[3,:]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )

U = array(zeros((4,N)))

U[0,0] = 1
U[1,0] = 0
U[2,0] = 0
U[3,0] = 1


def CN( dt, U, N ):

    for i in range(N-1):
        
        def CN_res(x):
            
            return x - U_temp - dt/2 * F_Kepler(x)
        U_temp = U[:,i] + dt/2 * F_Kepler(U[:,i])
        
        U[:,i+1] = fsolve(CN_res,U[:,i])
        
    return U
U = CN(dt, U, N)
plt.plot(U[0,N-1],U[1,N-1])
plt.show()