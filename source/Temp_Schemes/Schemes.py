from scipy.optimize import newton


def Euler(F,U,t0, tf):
    dt = tf-t0
    return U + dt * F(U,t0)

def RK4(F,U,t0, tf ):
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


