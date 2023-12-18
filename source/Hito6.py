from numpy import           array, zeros, linspace, sqrt, sin, cos, hstack, pi, meshgrid
import matplotlib.pyplot as plt
from Integrador import Integrate_ODE
from Temp_Schemes.Schemes import RK4
from ODEs.ODE import CR3BP
from scipy.optimize import newton

#--------------------------------------- CR3BP ---------------------------------------


t0 = 0     # Non dimensional initial time
tf = 20    # Non dimensional final time
#N = 10000
#t = linspace(0, N*dt, N+1)
dt = 0.01                       
t = linspace(t0, tf, int((tf-t0)/dt))
m1 = 5.972E24       # Mass of main body (Earth)
m2 = 7.348E22       # Mass of second body (Moon)
pi2 = m2 /(m1+m2)

# Initial conditions 
U0 = array(zeros((6,len(t)-1)))
U0 = array( [1-pi2, 0.0455, 0, -0.5, 0.5, 0] ) # x, y, z, vx, vy, vz

U = Integrate_ODE(U0, CR3BP, t, RK4)

#Plot trajectory in CR3BP
fig, ax = plt.subplots(figsize=(5,5), dpi=96)
# Moon orbit, plotted just to see magnitude, as in sinodical reference frame the moon orbit moves with the  x axis.
x_2 = (1 - pi2) * cos(linspace(0, pi, 100))
y_2 = (1 - pi2) * sin(linspace(0, pi, 100))

ax.plot( U[0,:], U[1,:], 'r', label ='Trajectory' )
ax.axhline(0, color='k')
ax.plot(hstack((x_2, x_2[::-1])), hstack((y_2, -y_2[::-1])))
ax.plot(-pi2, 0, 'bo', label="Earth")
ax.plot(1 - pi2, 0, 'go', label="Moon")
ax.plot(U0[0], U0[1], 'ro')
ax.legend(loc="upper left")
ax.set_aspect("equal")
plt.title("Arenstof orbit in CR3BP")
plt.show()

#------------------------------- Lagrange points -------------------------------

#Collinear points: L1, L2 and L3
def collinear_lagrange(x, pi_2):

    first_term = x
    second_term = (1 - pi_2) / abs(x + pi_2)**3 * (x + pi_2)
    third_term = pi_2 / abs(x - 1 + pi_2)**3 * (x - 1 + pi_2)
    return first_term - second_term - third_term

L_2 = newton(func=collinear_lagrange, x0=1, args=(pi2,))
L_1 = newton(func=collinear_lagrange, x0=0, args=(pi2,))
L_3 = newton(func=collinear_lagrange, x0=-1, args=(pi2,))
print(f"{L_1=}, {L_2=}, {L_3=}")

fig, ax = plt.subplots(figsize=(5,5), dpi=96)
ax.set_xlabel("$x^*$")
ax.set_ylabel("$y^*$")

ax.axhline(0, color='k')
ax.plot(hstack((x_2, x_2[::-1])), hstack((y_2, -y_2[::-1])))
ax.plot([-pi2, 0.5 - pi2, 1 - pi2, 0.5 - pi2, -pi2], [0, sqrt(3)/2, 0, -sqrt(3)/2, 0], 'k', ls="--", lw=1)

ax.plot(L_1, 0, 'rv', label="$L_1$") # L1
ax.plot(L_2, 0, 'r^', label="$L_2$") # L2
ax.plot(L_3, 0, 'rp', label="$L_3$") # L3
ax.plot(0.5 - pi2, sqrt(3)/2, 'rX', label="$L_4$") # L4
ax.plot(0.5 - pi2, -sqrt(3)/2, 'rs', label="$L_5$") # L5
ax.plot(-pi2, 0, 'bo', label="Earth")
ax.plot(1 - pi2, 0, 'go', label="Moon")
ax.legend(loc="upper left")
ax.set_aspect("equal")
plt.title("Lagrange points")
plt.show()


#-------------------------- Stability of Lagrange points ---------------------------------------

def Potential_function(x, y):
    phi = sqrt((x-1+pi2)**2+y**2)
    sigma = sqrt((x+pi2)**2+y**2)
    return -(1-pi2)/sigma - (pi2)/phi -0.5* ((1-pi2)*sigma**2 + pi2*phi**2)
	
x = linspace(-1.3, 1.3, 30)
y = linspace(-1.3,1.3, 30)

X, Y = meshgrid(x, y)
Z = Potential_function(X, Y)


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 100, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('3D contour')
plt.title("Stability contour of Lagrange points")
plt.show()


#-------------------------- Plot periodic orbits around Lagrange points --------------------------



t0 = 0     # Non dimensional initial time
tf = 1   # Non dimensional final time
dt = 0.0001 
t = linspace(t0, tf, int((tf-t0)/dt))

r12 = 384400 #km, Earth-Moon distance, used to erase dimension of coordinates
mu = 3.986e5
tc = sqrt(r12**3 / mu)

U0_Halo_L1 = array(zeros((6,len(t)-1)))
U0_Halo_L1 = array( [346485.4358330568/r12, 0, -2655.9162283251/r12, 0.02625013389512972*tc/r12, -0.01360917147008038*tc/r12, 0] )
U0_Lissajous_L1 = array( [345600.042099542/r12, 0, 0, 0.02611436483259207*tc/r12, -0.13252193918320*tc/r12, 0] )
U0_Halo_L2 = array( [476333.1026641181/r12, 0, -4478.74238999023/r12, 0.03701671556435333*tc/r12, -0.1541259694652506*tc/r12, 0] )
U0_Lissajous_L2 = array( [464183.4990431696/r12, 13912.28222140923/r12, 4519.357760825184/r12, 0.05782869445423985*tc/r12, 0.01817750924466563*tc/r12, 0.016320563401210*tc/r12] )

U_Halo_L1 = Integrate_ODE(U0_Halo_L1, CR3BP, t, RK4)
U_Lissajous_L1 = Integrate_ODE(U0_Lissajous_L1, CR3BP, t, RK4)

t0 = 0     # Non dimensional initial time
tf = 2   # Non dimensional final time
dt = 0.00001 
t = linspace(t0, tf, int((tf-t0)/dt))
""" U_Halo_L2 = Integrate_ODE(U0_Halo_L2, CR3BP, t, RK4)
U_Lissajous_L2 = Integrate_ODE(U0_Lissajous_L2, CR3BP, t, RK4) """
 
# L1 Halo and Lissajous orbits 2D
fig, ax = plt.subplots(figsize=(5,5), dpi=96)
ax.plot( U_Halo_L1[0,:], U_Halo_L1[1,:], 'r', label ='Halo L1' )
ax.plot( U_Lissajous_L1[0,:], U_Lissajous_L1[1,:], 'b', label ='Lissajous L1' )
ax.plot(hstack((x_2, x_2[::-1])), hstack((y_2, -y_2[::-1])))
ax.plot(-pi2, 0, 'bo', label="$Earth$")
ax.plot(1 - pi2, 0, 'go', label="Moon")
ax.legend(loc="upper left")
plt.title("Halo and Lissajous L1 orbits in sinodical frame")
plt.show()

# L2 Halo and Lissajous orbits 2D
""" fig, ax = plt.subplots(figsize=(5,5), dpi=96)
ax.plot( U_Halo_L2[0,:], U_Halo_L2[1,:], 'y', label ='Halo L2' )
ax.plot( U_Lissajous_L2[0,:], U_Lissajous_L2[1,:], 'g', label ='Lissajous L2' )
ax.plot(hstack((x_2, x_2[::-1])), hstack((y_2, -y_2[::-1])))
ax.plot(-pi2, 0, 'bo', label="$Earth$")
ax.plot(1 - pi2, 0, 'go', label="Moon")
ax.legend(loc="upper left")
ax.set_aspect("equal")
plt.show() """

# L1 Halo and Lissajous orbits 3D
ax = plt.axes(projection='3d')
ax.plot3D(U_Halo_L1[0,:], U_Halo_L1[1,:], U_Halo_L1[2,:], 'r', label ='Halo L1')
ax.plot3D(U_Lissajous_L1[0,:], U_Lissajous_L1[1,:], U_Lissajous_L1[2,:], 'b', label ='Lissajous L1')
ax.plot(1 - pi2, 0, 'go', label="Moon")
ax.legend(loc="upper left")
plt.title("3D Halo and Lissajous orbits")
plt.show()

# L2 Halo and Lissajous orbits 3D
""" ax = plt.axes(projection='3d')
ax.plot3D(U_Halo_L2[0,:], U_Halo_L2[1,:], U_Halo_L2[2,:], 'y', label ='Halo L2')
ax.plot3D(U_Lissajous_L2[0,:], U_Lissajous_L2[1,:], U_Lissajous_L2[2,:], 'g', label ='Lissajous L2')
ax.legend(loc="upper left")
plt.show() """