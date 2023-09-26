from numpy import array, zeros
import matplotlib.pyplot as plt

def F_Kepler(U):
    
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )

def F_Kepler_2():
N = 10000
U = array( [1, 0, 0, 1] )
dt = 0.001
x = array( zeros(N) )
y = array( zeros(N) )
vx = array( zeros(N) )
vy = array( zeros(N) )

x[0] = U[0]
y[0] = U[1]
vx[0] = U[2]
vy[0] = U[3]

for i in range(0, N):
    F = F_Kepler( U )
    U = U + dt/2 * ( F )
    x[i] = U[0]
    y[i] = U[1]
    
plt.plot(x, y)
plt.show()

import numpy
import matplotlib.pyplot as plt

# Discretización en el tiempo
num_steps = 1000
dt = 0.01
times = np.linspace(0, T, num_steps)

# Inicialización de arrays para almacenar los resultados
x = zeros(num_steps)
y = zeros(num_steps)
vx = zeros(num_steps)
vy = zeros(num_steps)

# Condiciones iniciales
x[0] = 1
y[0] = 0.0
vx[0] = 0.0
vy[0] = 1

# Método de Crank-Nicolson para la integración numérica
for i in range(1, num_steps):
    r = sqrt(x[i-1]**2 + y[i-1]**2)
    ax_n = - (G * M * x[i-1]) / (r**3)
    ay_n = - (G * M * y[i-1]) / (r**3)

    # Aproximación de la derivada temporal usando Crank-Nicolson
    ax_n1 = - (G * M * x[i]) / (sqrt(x[i]**2 + y[i]**2)**3)
    ay_n1 = - (G * M * y[i]) / (sqrt(x[i]**2 + y[i]**2)**3)

    # Actualización de posiciones y velocidades usando Crank-Nicolson
    x[i] = x[i-1] + dt * 0.5 * (vx[i-1] + vx[i] + dt * 0.5 * (ax_n + ax_n1))
    y[i] = y[i-1] + dt * 0.5 * (vy[i-1] + vy[i] + dt * 0.5 * (ay_n + ay_n1))
    vx[i] = vx[i-1] + dt * 0.5 * (ax_n + ax_n1)
    vy[i] = vy[i-1] + dt * 0.5 * (ay_n + ay_n1)

# Graficar la órbita
plt.figure(figsize=(8, 6))
plt.plot(x, y)
plt.xlabel('Coordenada x')
plt.ylabel('Coordenada y')
plt.title('Órbita Kepleriana en 2D (método de Crank-Nicolson)')
plt.axis('equal')
plt.grid(True)
plt.show()





