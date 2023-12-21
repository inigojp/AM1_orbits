
from numpy import array, zeros, log10, ones, vstack  
from numpy.linalg import norm, lstsq
from Temp_Schemes.Schemes import *
from ODEs.ODE import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton
from Error.Richardson import *
from Integrador import Integrate_ODE

#---------------------------------------- Richardson, calculo de error ----------------------------------------


#Inputs
t0 = 0
tf = 50
dt = 0.01

# Discretización de tiempo
N = int((tf-t0) / dt)                    
t = linspace(t0, tf, N)  # Tiempo inicial, tiempo final, número de puntos

# Initial conditions 

U0_c = array( [1, 0] )

U0_k = array( [1, 0, 0, 1] )

order_RK4 = 4
order_CN = 2
order_E = 1
order_IE = 1

Richardson_output_RK4 = Error_ODE( t, order_RK4, U0_c, ODE_Rayleigh, RK4 )
Richardson_output_CN = Error_ODE( t, order_CN, U0_c, ODE_Rayleigh, Crank_Nicolson )
Richardson_output_IE = Error_ODE( t, order_IE, U0_c, ODE_Rayleigh, Inverse_Euler ) 
Richardson_output_E = Error_ODE( t, order_E, U0_c, ODE_Rayleigh, Euler )

Error_RK4 = Richardson_output_RK4[0]
Error_CN = Richardson_output_CN[0]
Error_IE = Richardson_output_IE[0]
Error_E = Richardson_output_E[0]

#Plot error para cada instante de integracion para la coordenada x con difererntes metodos
plt.plot(log10(Error_RK4[0,:]), label = "RK4")
plt.plot(log10(Error_CN[0,:]), label = "Crank Nicolson")
plt.plot(log10(Error_IE[0,:]), label = "Inverse Euler")
plt.plot(log10(Error_E[0,:]), label = "Euler")
plt.legend(loc="upper left")
plt.title("Error temporal obtenido por extrapolación de Richardson")
plt.show()


#----------------------------------- Convergencia ----------------------------------------
N_mallas = 2 # El tiempo de computo es extremedamente dependiente de el número de mallas, pueidendo alargarse a minutos para N_mallas > 5, tf > 50 y dt < 0.01
Convergence_RK4 = Temporal_convergence_rate( t, N_mallas, U0_k, F_Kepler, RK4 )
Convergence_CN = Temporal_convergence_rate( t, N_mallas, U0_k, F_Kepler, Crank_Nicolson )
""" Convergence_IE = Temporal_convergence_rate( t, N_mallas, U0_k, F_Kepler, Inverse_Euler ) """
Convergence_E = Temporal_convergence_rate( t, N_mallas, U0_k, F_Kepler, Euler )

log_E_RK4 = abs(Convergence_RK4[1])
log_N_RK4 = Convergence_RK4[2]

log_E_CN = abs(Convergence_CN[1])
log_N_CN = Convergence_CN[2]

""" log_E_IE = abs(Convergence_IE[1])
log_N_IE = Convergence_IE[2]  """

log_E_E = abs(Convergence_E[1])
log_N_E = Convergence_E[2] 

# Plot de error de error para diferentes mallas, a mayor numero de malla menor time step (malla mas refinada)
plt.plot(log_N_RK4[:], log_E_RK4[:], label = "RK4")
plt.plot(log_N_CN[:], log_E_CN[:], label = "Crank-Nicolson")
""" plt.plot(log_N_IE[:], log_E_IE[:], label = "Inverse-Euler") """
plt.plot(log_N_E[:], log_E_E[:], label = "Euler")
plt.legend(loc="upper left")
plt.show()

