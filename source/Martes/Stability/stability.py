from numpy import array, zeros, linspace, abs, transpose, float64
import matplotlib.pyplot as plt

def Euler(F,i):
    
    return U[:,i] + dt * F(U[:,i])
    
def Crank_Nicolson(F,i):
    
    def g(x):
        return x - a -dt/2 * (F(a) + F(x))
    a = U[:,i]
    
    return newton(g, U[:,i])

def RK4(U,,dt,F,i):
    
    k1 = F( U[:,i], t )
    k2 = F( U[:,i] +k1*dt/2, t + dt/2 )
    k3 = F( U[:,i] +k2*dt/2, t + dt/2 )
    k4 = F( U[:,i] + k3*dt, t + dt )
     
    k = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    return U[:,i] + dt * k

def Inverse_Euler(F,i):
    def g(x):
        return x - a - dt * F(x,t) 
    a = U[:,i]
    
    return newton(g, U[:,i])


def Stability_Region(Scheme, N, x0, xf, y0, yf):

    x, y = linspace(x0, xf, N), linspace(y0, yf, N)
    rho = zeros( (N,N), dtype=float64)

    for i in range(N):
        for j in range(N):

            w = complex(x[i], y[i])
            r = Scheme(1., 1., lambda u, t: w*u )
            rho[i, j] = abs(r)

        return rho, x, y
    
def test_Stability_regions():

    schemes = [RK4]

    for scheme in schemes:
        rho, x, y = Stability_Region(scheme, 100, -4, 2, -4, 4)
        plt.contour( x, y, transpose(rho), linspace(0, 1, 11) )
        plt.axis('equal')
        plt.grid()
        plt.show()
    
if __name__ == '__main__': 
    test_Stability_regions()

# if __name__ == '__main__': Sirve para que este codigo solo se ejecute si es el main, 
# de manera que si se importa este archivo en otro, al correr el otro este no se ejecute y muestre las gráficas