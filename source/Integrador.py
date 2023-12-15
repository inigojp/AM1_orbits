from numpy import array, zeros, log10, ones, vstack  
from numpy.linalg import norm, lstsq
from time import process_time


def Integrate_ODE(U0, F, t, scheme):
    N, Nv=  len(t)-1, len(U0)           # N = número evaluaciones, Nv = Número variables (x,y,vx,vy...)
    U = zeros( (Nv, N))
    U[:,0] = U0

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

def Cauchy_problem( F, t, U0, Temporal_scheme): 


 start = process_time()
 N, Nv=  len(t)-1, len(U0)
 U = zeros( (N+1, Nv), dtype=float) 

 U[0,:] = U0
 for n in range(N):

     U[n+1,:] = Temporal_scheme( U[n, :], t[n+1] - t[n], t[n],  F ) 

 finish = process_time()
 print("Cauchy_Problem, CPU Time=",finish-start," seconds.") 
 return U