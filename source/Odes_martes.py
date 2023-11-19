
from numpy import array, zeros, linspace
from Temp_Schemes.Schemes import *
import matplotlib.pyplot as plt
from ODEs.ODE import *

## Integración sistemas dinámicos

##Integrar 
def Integrate_ODE(U, F, t, scheme):

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

N = 30000
d_t = 0.01                        
t = linspace(0,N*d_t, N+1)


U2 = array(zeros((2,len(t)-1)))
U2[:,0] = array( [1, 0.01] )

U2 = Integrate_ODE(U2, ODE_Rayleigh2, t, Inverse_Euler)


f1 = plt.figure()
plt.plot(U2[0,10000:N],U2[1,10000:N])
f2 = plt.figure()
plt.plot(U2[0,0000:N])
#plt.plot(U2[1,5000:N])



plt.show()