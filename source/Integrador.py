from numpy import array, zeros, log10, ones, vstack  
from numpy.linalg import norm, lstsq
from time import process_time


def Integrate_ODE(U0, F, t, scheme):
    N = len(t)-1; Nv = len(U0)           # N = número evaluaciones, Nv = Número variables (x,y,vx,vy...)
    U = zeros( (Nv, N+1))
    U[:,0] = U0

    for i in range(N):
        U[:,i+1] = scheme(F, U[:,i], t[i+1] - t[i], t[i])
        
    return U

