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

def ODE_Van_Der_Por_2(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=-5*(x**2-1)*y-x

    return array( [vx, vy] )

def ODE_Rayleigh(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=x-x**3-y+0.6*math.cos(t)

    return array( [vx, vy] )

def ODE_Rayleigh2(U, t): 
    x, y = [U[0], U[1]]

    vx=y
    vy=x-x**3-y+3*math.cos(t)

    return array( [vx, vy] )

