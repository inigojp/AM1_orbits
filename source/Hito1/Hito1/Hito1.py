from numpy import array, zeros
import matplotlib.pyplot as plt

def F_Kepler(U):
    
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )

N = 10000
U = array( [1, 0, 0, 1] )
dt = 0.001
x = array( zeros(N) )
y = array( zeros(N) )
x[0] = U[0]
y[0] = U[1]

for i in range(0, N):
    F = F_Kepler( U )
    U = U + dt * F
    x[i] = U[0]
    y[i] = U[1]
    
plt.plot(x, y)
plt.show()
