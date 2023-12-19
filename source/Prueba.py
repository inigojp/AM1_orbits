from numpy import array, zeros, log10, ones, vstack  
from numpy.linalg import norm, lstsq
from Temp_Schemes.Schemes import *
from ODEs.ODE import *
import matplotlib.pyplot as plt
from time import process_time
from matplotlib.animation import FuncAnimation


def Cauchy_problem( F, t, U0, Temporal_scheme): 


 start = process_time()
 N, Nv=  len(t)-1, len(U0)
 U = zeros( (N+1, Nv), dtype=float) 

 U[0,:] = U0
 for n in range(N):

     U[n+1,:] = Temporal_scheme(F, U[n, :], t[n], t[n+1] ) 

 finish = process_time()
 print("Cauchy_Problem, CPU Time=",finish-start," seconds.") 
 return U

N = 18000
d_t = 0.01                        
t = linspace(0,N*d_t, N+1)

U0 = array(zeros((len(t)-1,4)))
U0 = array( [0.994, 0, 0, -2.001585106] )
#U0 = array(zeros((len(t)-1,4)))
#U0 = array( [1, 0, 0, 1] )

U = Cauchy_problem(Arenstorf, t, U0, leap_frog)


plt.plot(U[:,0], U[:,1])
plt.show()