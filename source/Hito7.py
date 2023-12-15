from numpy import array, zeros, linspace, arange, min
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton
from Temp_Schemes.Schemes import *
from ODEs.ODE import Arenstorf
from matplotlib.animation import FuncAnimation

def Integrate_ODE(U, F, t, scheme):

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

N = 18000
d_t = 0.001                        
t = linspace(0,N*d_t, N+1)


U0 = array(zeros((4,len(t)-1)))
U0[:,0] = array( [0.994, 0, 0, -2.001585106] )


U = Integrate_ODE(U0, Arenstorf, t, GBS_Scheme)


fig, ax = plt.subplots()
ax.set_xlim(min(U[0]) - 1, max(U[0]) + 1)
ax.set_ylim(min(U[1]) - 1, max(U[1]) + 1)
line, = ax.plot([], [], 'b-')
points, = ax.plot([], [], 'bo')  # Puntos azules

# Función de inicialización
def init():
    line.set_data([], [])
    points.set_data([], [])
    return line,points

# Función de animación
def update(frame):
    x = U[0, :frame+1]
    y = U[1, :frame+1]
    line.set_data(x, y)
    points.set_data(U[0, frame], U[1, frame])
    return line, points

# Crear la animación
num_frames = U.shape[1] 
ani = FuncAnimation(fig, update, frames=num_frames, init_func=init, interval=0.0001, blit=True)

plt.show()