from numpy import array
import math

def F_Cauchy(U, t=0): # df/dx = x

    x = U[0]
    return x

def ODE(U, t=0): # df/dx = x
    x, y = [U[0], U[1]]
    vx=1-4*x+x**2*y**2
    vy=3*x-x**2*y

    return array( [vx, vy] )

def F_Kepler(U, t=0):
 
    x, y, vx, vy = [U[0], U[1], U[2], U[3]]
    mr = (x**2 + y**2) ** 1.5
    
    return array( [vx, vy, -x/mr, -y/mr] )

def ODE_DUffing(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=-x**3-0.01*y+math.cos(t)

    return array( [vx, vy] )

def ODE_Van_Der_Por_1(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=-0.2*(x**2-1)*y-x

    return array( [vx, vy] )

def ODE_Rayleigh(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=x-x**3-y+0.6*math.cos(t)

    return array( [vx, vy] )

def Lorentz(U, t):
    x, y, z =  [U[0], U[1], U[2]]

    vx = -10*x+10*y
    vy = 28*x-x*z-y
    vz = -(8/3)*z+x*y

    return array( [vx, vy, vz] )

def Rossler(U, t):
    x, y, z =  [U[0], U[1], U[2]]
    #a, b, c = 0.398,2, 4
    vx = -y-z
    vy = 0.398*y+x
    vz = 2+z*(x-4)

    return array( [vx, vy, vz] )

def Rossler2(U, t):
    x, y, z =  [U[0], U[1], U[2]]
    #a, b, c = 0.398,2, 4
    vx = -y-z
    vy = 0.17*y+x
    vz = 0.4+z*(x-8.5)

    return array( [vx, vy, vz] )



