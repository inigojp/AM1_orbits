from numpy import array, linspace, min
import matplotlib.pyplot as plt
from Temp_Schemes.Schemes import *
from ODEs.ODE import Arenstorf
from matplotlib.animation import FuncAnimation
from Integrador import *

# Inputs
t0 = 0
tf = 17.065216 
dt = 0.01

# Discretización de tiempo
N = int((tf-t0) / dt)                    
t = linspace(t0, tf, N)  # Tiempo inicial, tiempo final, número de puntos
NL_fixed = 1

# Arenstof initial conditions
U0 = array( [0.994, 0, 0, -2.001585106] )

U = Integrate_ODE_GBS(U0,  Arenstorf, t, GBS_Scheme, NL_fixed)

plt.plot(U[:,0], U[:,1])
plt.show()


# Creacion de gif
""" anim_U = U
fps = 30
duration = 5
anim_save_name = "Arenstorf_anim.gif"

def create_animation(anim_U, anim_fps, anim_duration, anim_save_name,N):
   # Funciones para crear la animacion:
   # - Función de inicializacion: se llama para crear la trama base vacia
   def init():
      line.set_data([], [])
      return line,
   # - Funcion de actualizacion para cada cuadro de la animacion
   def update(frame):
      x1 = anim_U[0, :frame]  # Obtener datos hasta el cuadro actual
      y1 = anim_U[1, :frame]
      line.set_data(x1, y1)
      return line,
 
   # Crear la figura usando plt.figure
   plt.figure(figsize=(8, 8))
   plt.tick_params(labelsize=13)
 
   # Crear una linea vacia que se actualizara en la animacion
   line, =  plt.plot([], [], '-',color='blue')
 
   # Definir los parámetros de la animación
   anim_frames = linspace(0, N, int(anim_fps*anim_duration+1), dtype=int)
   anim_interval = (anim_duration / int(anim_fps*anim_duration+1)) * 1000
 
   # Definir los límites de los ejes de la grafica
   rangex = (max(anim_U[0,:]) - min(anim_U[0,:])) * 0.1
   rangey = (max(anim_U[1,:]) - min(anim_U[1,:])) * 0.1
   xlim = [min(anim_U[0,:]) - rangex, max(anim_U[0,:]) + rangex]
   ylim = [min(anim_U[1,:]) - rangey, max(anim_U[1,:]) + rangey]
 
   # Definicion de propiedades de la grafica
   plt.xlim(xlim)
   plt.ylim(ylim)
   plt.xlabel('$x$', fontdict = {'fontsize':14, 'fontweight':'normal', 'color':'k'})
   plt.ylabel('$dx/dt$', fontdict = {'fontsize':14, 'fontweight':'normal', 'color':'k'})
   plt.plot(0.1,0,marker='o', color='b')
 
   # Creacion y guardado de la animacion
   ani = FuncAnimation(plt.gcf(), update, frames=anim_frames, init_func=init, blit=True, interval=anim_interval)
   ani.save(anim_save_name, writer='Pillow', fps=anim_fps)
   plt.show()
Animation = create_animation(U, fps, duration, anim_save_name, N)
 """