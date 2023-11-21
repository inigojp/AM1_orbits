
from numpy import array, zeros, linspace
from Temp_Schemes.Schemes import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from ODEs.ODE import *

## Integración sistemas dinámicos

##Integrar 
def Integrate_ODE(U, F, t, scheme):

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

N = 10000
d_t = 0.01                        
t = linspace(0,N*d_t, N+1)


U0 = array(zeros((3,len(t)-1)))
U0[:,0] = array( [0.5, 1, 2] )
U1 = array(zeros((3,len(t)-1)))
U1[:,0] = array( [0.5, 1, 2] )

U = Integrate_ODE(U0, Rossler, t, RK4)
U2 = Integrate_ODE(U1, Rossler2, t, RK4)

f1 = plt.figure()
ax = plt.axes(projection ='3d')
ax.plot3D(U[0,:],U[1,:],U[2,:])
f2 = plt.figure()
ax = plt.axes(projection ='3d')
ax.plot3D(U2[0,:],U2[1,:],U2[2,:])
#plt.plot(U2[1,5000:N])



plt.show()