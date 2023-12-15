from scipy.optimize import newton
from numpy import array, zeros, linspace, arange, min, dot

def Euler(F, U, t0, tf):
    dt = tf-t0
    return U + dt * F(U,t0)

def RK4(F, U, t0, tf ):
    dt = tf-t0
    k1 = F( U, t0 )
    k2 = F( U +k1*dt/2, t0 + dt/2 )
    k3 = F( U +k2*dt/2, t0 + dt/2 )
    k4 = F( U + k3*dt, t0 + dt )
     
    k = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    return U + dt * k


def Crank_Nicolson(F, U, t0, tf):
    dt = tf-t0
    def g(x):
        return x - a -dt/2 * (F(a,t0) + F(x,tf))
    a = U
    
    return newton(g, U)

def Inverse_Euler(F, U, t0, tf):
    dt = tf-t0
    def g(x):
        return x - a - dt * F(x,tf) 
    a = U
    
    return newton(g, U)

#Schemes adapted for stability 
def Euler_st( U, dt, t, F): 

   return U + dt * F( U, t ) 

def RK4_st( U, dt, t, F): 

     k1 = F( U, t)
     k2 = F( U + dt * k1/2, t + dt/2 )
     k3 = F( U + dt * k2/2, t + dt/2 )
     k4 = F( U + dt * k3,   t + dt   )
 
     return U + dt * ( k1 + 2*k2 + 2*k3 + k4 )/6

def Inverse_Euler_st( U, dt, t, F): 

   def  Residual_IE(X):  

          return  X - U - dt * F(X, dt + t) 

   return  newton( Residual_IE, U )
    

def Crank_Nicolson_st(U, dt, t, F ): 

    def Residual_CN(X): 
         
         return  X - a - dt/2 *  F(X, t + dt)

    a = U  +  dt/2 * F( U, t)  
    return newton( Residual_CN, U )

#  GBS
Tolerance = 1e-8
NL_fixed = 3
N_GBS_effort = 0

def GBS_Scheme(F, U1, t0, tf):
    global Tolerance
    global NL_fixed
    global N_GBS_effort

    if NL_fixed > 0:
        GBS_solution_NL(U1, t0, tf, F, NL_fixed)
    else:
        raise ValueError("NL_fixed must be greater than 0 for fixed scheme")

def GBS_solution_NL(U1, t0, tf, F, NL):
    N = mesh_refinement(NL)
    U = zeros((NL, len(U1)))

    for i in range(NL):
        Modified_midpoint_scheme(U1, t0, tf, F, U[i, :], N[i])

    U2 = Corrected_solution_Richardson(N, U)
    return U2

def mesh_refinement(Levels):
    return arange(1, Levels + 1)

def Modified_midpoint_scheme(U1, t0, tf, F, U2, N):
    h = (tf - t0) / (2 * N)
    U = zeros((len(U1), 2 * N + 2))
    U[:, 0] = U1
    U[:, 1] = U[:, 0] + h * F(U[:, 0], t0)

    for i in range(1, 2 * N + 1):
        ti = t0 + i * h
        U[:, i + 1] = U[:, i - 1] + 2 * h * F(U[:, i], ti)

    U2[:] = (U[:, 2 * N + 1] + 2 * U[:, 2 * N] + U[:, 2 * N - 1]) / 4.0

def Corrected_solution_Richardson(N, U):
    h = 1.0 / (2 * N)
    x = h ** 2
    W = 1.0 / (x * (x - 1))
    Uc = dot(W / x, U) / sum(W / x)
    return Uc
