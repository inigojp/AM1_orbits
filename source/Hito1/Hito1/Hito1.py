from numpy import array, zeros, linspace
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

N = 1000                         # PASOS DE INTEGRACION
dt = 0.01                        # INTERVALO ENTRE PASOS
t_values = array(zeros(N))

def F_Kepler(t, U):
 
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )

#EULER
Ue = array( [1, 0, 0, 1] )
x = array( zeros(N) )
y = array( zeros(N) )
x[0] = Ue[0]
y[0] = Ue[1]

for i in range(0, N):
    t = t_values[i]
    F = F_Kepler(t, Ue )
    Ue = Ue + dt * F
    x[i] = Ue[0]
    y[i] = Ue[1]
    
plt.plot(x, y)
plt.show()

# CRANK-NICOLSON

# U_cn  = array([1.0, 0.0, 0.0, 1.0])
# dt_cn = 0.01

# x_cn = array(zeros(N))
# y_cn = array(zeros(N))
# x_cn[0] = U_cn[0]
# y_cn[0] = U_cn[1]

# for i in range(1, N):
#     U_n = U_cn                              # Valor en el tiempo tn
#     F_n = F_Kepler(t,U_n)                     # Derivadas en el tiempo tn
    
#     U_half = U_n + 0.5 * dt_cn * F_n        # Valor en el tiempo tn+1/2
#     F_half = F_Kepler(t,U_half)               # Derivadas en el tiempo tn+1/2
    
#     U_cn = U_n + dt_cn * F_half             # Valor en el tiempo tn+1
#     x_cn[i] = U_n[0]
#     y_cn[i] = U_n[1]
    
# plt.plot(x_cn, y_cn)
# plt.show()
U = array(zeros((4,N)))
x = array( zeros(N) )
y = array( zeros(N) )
U[0,0] = 1
U[1,0] = 0
U[2,0] = 0
U[3,0] = 1
x[0] = U[0,0]
y[0] = U[1,0]

def F_Kepler(U):
 
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )

def CN( dt, U, N ):

    for i in range(N-1):
        
        def CN_res(x):
            
            return x - U_temp - dt/2 * F_Kepler(x)
        U_temp = U[:,i] + dt/2 * F_Kepler(U[:,i])
        
        U[:,i+1] = fsolve(CN_res,U[:,i])
        x[i+1] = U[0,i+1]
        y[i+1] = U[1,i+1]
    return U
U = CN(dt, U, N)
plt.plot(x,y)
plt.show()

#RUNGE-KUTTA 4TH ORDER
def F_Kepler(t, U):
 
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )
U_rk  = array([1.0, 0.0, 0.0, 1.0])

x_rk = array(zeros(N))
y_rk = array(zeros(N))
x_rk[0] = U_rk[0]
y_rk[0] = U_rk[1]

for i in range(0,N):
    t = t_values[i]
    x_rk[i] = U_rk[0]
    y_rk[i] = U_rk[1]
    k1 = F_Kepler( t, U_rk )
    k2 = F_Kepler( t + dt/2, U_rk +k1*dt/2 )
    k3 = F_Kepler( t + dt/2, U_rk +k2*dt/2 )
    k4 = F_Kepler( t + dt, U_rk + k3*dt )
     
    k = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    U_rk = U_rk + dt * k
    
plt.plot(x_rk, y_rk)
plt.show()