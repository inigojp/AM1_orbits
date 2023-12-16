from numpy import array, zeros, log10, ones, vstack  
from numpy.linalg import norm, lstsq
from Temp_Schemes.Schemes import *
from ODEs.ODE import F_Kepler
import matplotlib.pyplot as plt

def Integrate_ODE(U0, F, t, scheme):
    N, Nv=  len(t)-1, len(U0)
    U = zeros( (Nv, N), dtype=float)
    U[:,0] = U0

    for i in range(0, N-1):
        U[:,i+1] = scheme(F,U[:,i],t[i],t[i+1])
        
    return U

t0 = 0
tf = 1000
dt = 0.01         
Np = int((tf-t0)/dt)               
t = linspace(t0, tf, Np)


U0_k = array( [1, 0, 0, 1] )

U_k = Integrate_ODE(U0_k, F_Kepler, t, RK4)

plt.plot(U_k[0,:],U_k[1,:])
plt.axis('equal')
plt.show()




r12 = 384400 #km, Earth-Moon distance, used to erase dimension of coordinates
mu = 3.986e5
tc = sqrt(r12**3 / mu)

U0_Halo_L1 = array(zeros((6,len(t)-1)))
U0_Halo_L1 = array( [346485.4358330568/r12, 0, -2655.9162283251/r12, 0.02625013389512972*tc/r12, -0.01360917147008038*tc/r12, 0] )

U_Halo_L1 = Integrate_ODE(U0_Halo_L1, CR3BP, t, RK4)

fig, ax = plt.subplots(figsize=(5,5), dpi=96)
ax.plot( U_Halo_L1[0,:], U_Halo_L1[1,:], 'r', label ='Trajectory' )
plt.show()


def Embedded_RK(F, U, t0, tf, t, q, Tolerance): 
    dt = tf-t0
    #(a, b, bs, c) = Butcher_array(q)
    #a, b, bs, c = Butcher_array(q)
 
    N_stages = { 2:2, 3:4, 8:13  }
    Ns = N_stages[q]
    a = zeros( (Ns, Ns), dtype=float) 
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


def RK_stages( F, U, t0, tf, t, a, c ): 
     dt = tf-t0
     k = zeros( (len(c), len(U)), dtype=float )

     for i in range(len(c)): 

        for  j in range(len(c)-1): 
          Up = U + dt * dot( a[i, :], k)

        k[i, :] = F( Up, t + c[i] * dt ) 

     return k 