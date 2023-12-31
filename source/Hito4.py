from numpy import array, zeros, log10, ones, vstack, transpose
from numpy.linalg import norm, lstsq
from Temp_Schemes.Schemes import *
from ODEs.ODE import Oscillator
from Integrador import Integrate_ODE
from scipy.optimize import fsolve, newton
from Error.Stability import Stability_Region
import matplotlib.pyplot as plt

#----------------- Oscilador ----------------------------------
#Inputs
t0 = 0
tf = 50
dt = 0.01

# Discretización de tiempo
N = int((tf-t0) / dt)                    
t = linspace(t0, tf, N)  # Tiempo inicial, tiempo final, número de puntos

#Initial conditions Oscillator

schemes = [RK4, Crank_Nicolson, Euler, Inverse_Euler] 
U0_o = array( [1, 5] )
for scheme in schemes: 
    U_o = Integrate_ODE(U0_o, Oscillator, t, scheme)
    plt.plot(U_o[0,:],U_o[1,:])
    plt.axis('equal')
    plt.grid()
    plt.show()

#------------------ Region stabilidad numerica -----------------------

schemes = [RK4, Crank_Nicolson, Euler, Inverse_Euler] 
Ns = 100
for scheme in schemes: 
    rho, x, y  = Stability_Region(scheme, Ns, -4, 2, -4, 4)
    plt.contour( x, y, transpose(rho), linspace(0, 1, 11) )
    plt.axis('equal')
    plt.grid()
    plt.show()