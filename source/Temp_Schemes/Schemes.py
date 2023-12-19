from matplotlib.pylab import eig
from scipy.optimize import newton
from numpy import array, eye, zeros, linspace, arange, min, dot, float64
from numpy.linalg import norm

def Euler(F, U, dt, t ):

    return U + dt * F(U,t)

def RK4(F, U, dt, t ):

    k1 = F( U, t )
    k2 = F( U +k1*dt/2, t + dt/2 )
    k3 = F( U +k2*dt/2, t + dt/2 )
    k4 = F( U + k3*dt, t + dt )
     
    k = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    return U + dt * k


def Crank_Nicolson(F, U, dt, t):
 
    def g(x):
         return  x - a - dt/2 *  F(x, t + dt)

    a = U  +  dt/2 * F( U, t) 
    
    return newton(g, U)

def Inverse_Euler(F, U, dt, t):

    def g(x):
        return x - a - dt * F(x,t) 
    a = U
    
    return newton(g, U)


# ------------------------------------------------- GBS -----------------------------------------------
Tolerance = 1e-8
NL_fixed = 3
N_GBS_effort = 1

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

#------------------------------------------ RK4 de paso variable ------------------------------------
def Embedded_RK( U, dt, t, F, q, Tolerance): 
    
    #(a, b, bs, c) = Butcher_array(q)
    #a, b, bs, c = Butcher_array(q)
 
    N_stages = { 2:2, 3:4, 8:13  }
    Ns = N_stages[q]
    a = zeros( (Ns, Ns), dtype=float64) 
    b = zeros(Ns); bs = zeros(Ns); c = zeros(Ns) 
   
    if Ns==2: 
     
     a[0,:] = [ 0, 0 ]
     a[1,:] = [ 1, 0 ] 
     b[:]  = [ 1/2, 1/2 ] 
     bs[:] = [ 1, 0 ] 
     c[:]  = [ 0, 1]  

    elif Ns==13: 
       c[:] = [ 0., 2./27, 1./9, 1./6, 5./12, 1./2, 5./6, 1./6, 2./3 , 1./3,   1., 0., 1.]

       a[0,:]  = [ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0] 
       a[1,:]  = [ 2./27, 0., 0., 0., 0., 0., 0.,  0., 0., 0., 0., 0., 0] 
       a[2,:]  = [ 1./36 , 1./12, 0., 0., 0., 0., 0.,  0.,0., 0., 0., 0., 0] 
       a[3,:]  = [ 1./24 , 0., 1./8 , 0., 0., 0., 0., 0., 0., 0., 0., 0., 0] 
       a[4,:]  = [ 5./12, 0., -25./16, 25./16., 0., 0., 0., 0., 0., 0., 0., 0., 0]
       a[5,:]  = [ 1./20, 0., 0., 1./4, 1./5, 0., 0.,0., 0., 0., 0., 0., 0] 
       a[6,:]  = [-25./108, 0., 0., 125./108, -65./27, 125./54, 0., 0., 0., 0., 0., 0., 0] 
       a[7,:]  = [ 31./300, 0., 0., 0., 61./225, -2./9, 13./900, 0., 0., 0., 0., 0., 0] 
       a[8,:]  = [ 2., 0., 0., -53./6, 704./45, -107./9, 67./90, 3., 0., 0., 0., 0., 0] 
       a[9,:]  = [-91./108, 0., 0., 23./108, -976./135, 311./54, -19./60, 17./6, -1./12, 0., 0., 0., 0] 
       a[10,:] = [ 2383./4100, 0., 0., -341./164, 4496./1025, -301./82, 2133./4100, 45./82, 45./164, 18./41, 0., 0., 0] 
       a[11,:] = [ 3./205, 0., 0., 0., 0., -6./41, -3./205, -3./41, 3./41, 6./41, 0., 0., 0]
       a[12,:] = [ -1777./4100, 0., 0., -341./164, 4496./1025, -289./82, 2193./4100, 51./82, 33./164, 19./41, 0.,  1., 0]
      
       b[:]  = [ 41./840, 0., 0., 0., 0., 34./105, 9./35, 9./35, 9./280, 9./280, 41./840, 0., 0.] 
       bs[:] = [ 0., 0., 0., 0., 0., 34./105, 9./35, 9./35, 9./280, 9./280, 0., 41./840, 41./840]     
     

    
    k = RK_stages( F, U, t, dt, a, c ) 
    Error = dot( b-bs, k )

    dt_min = min( dt, dt * ( Tolerance / norm(Error) ) **(1/q) )
    N = int( dt/dt_min  ) + 1
    h = dt / N
    Uh = U.copy()

    for i in range(0, N): 

        k = RK_stages( F, Uh, t + h*i, h, a, c ) 
        Uh += h * dot( b, k )

    return Uh


def RK_stages( F, U, t, dt, a, c ): 

     k = zeros( (len(c), len(U)), dtype=float64 )

     for i in range(len(c)): 

        for  j in range(len(c)-1): 
          Up = U + dt * dot( a[i, :], k)

        k[i, :] = F( Up, t + c[i] * dt ) 

     return k 
# ======================================== GBS ======================================
def leap_frog(F, U, t1, t2):
    dt = t2 - t1
    N = len(U)
    t_old = 0
    
    if t1 < t_old or t_old == 0:
        U2 = U + dt * F(U, t1)
    else:
        U2 = U0 + 2 * dt * F(U, t1)
        U0 = U
    
    t_old = t2

    return U2

def GBS_Scheme(U1, t1, t2, F, NL_fixed):
    N_steps = 0
    U1s = U1.copy()
    t1s = t1
    t2s = t1

    while t2s < t2:
        lambda_vals = Eigenvalues_Jacobian(F, U1s, t1s)
        dt = 0.1 / max(abs(lambda_vals))

        if t1s + dt > t2:
            t2s = t2
        else:
            t2s = t1s + dt

        U1s = GBS_Solution(F, t1s, t2s, U1s, NL_fixed)
        t1s = t2s
        N_steps += 1

    return U1s

def Eigenvalues_Jacobian(F, U, t, epsilon=1e-8):
    N = len(U)
    identity_matrix = eye(N)
    jacobian_matrix = zeros((N, N))

    # Calcula las derivadas parciales numericamente
    for i in range(N):
        U_perturbed = U.copy()
        U_perturbed[i] += epsilon
        F_perturbed = F(U_perturbed, t)
        jacobian_matrix[:, i] = (F_perturbed - F(U, t)) / epsilon

    # Usa eig de scipy para obtener los valores propios
    eigenvalues, _ = eig(jacobian_matrix)

    return eigenvalues

def GBS_Solution(F, t1, t2, U1, Tolerance=1e-8):
    NL = 1
    Error = 1.0
    NLmax = 8

    Ucs = U1.copy()
    Uc = U1.copy()
    UL = zeros((NLmax+1, len(U1)))

    while Error > Tolerance and NL < NLmax:
        NL += 1
        GBS_solutionL(F, t1, t2, U1, UL, Uc, Ucs, NL)
        Error = norm(Uc - Ucs)

    return Ucs

def GBS_solutionL(F, t1, t2, U1, UL, Uc, Ucs, NL):
    lambda_vals = Eigenvalues_Jacobian(F, U1, t1)
    dt = (t2 - t1) / NL
    next_level = True

    for i in range(1, NL + 1):
        UL[i, :] = leap_frog(U1, t1, t1 + i * dt, F)

        if not next_level:
            continue

        if i > 1:
            Uc[:] = Ucs + (UL[i, :] - Ucs) / (2 ** (2 * i) - 1)

            if norm(Uc - Ucs) > 0.00001:
                next_level = False

    Ucs[:] = Uc