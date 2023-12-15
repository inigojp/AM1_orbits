from numpy import array, zeros, log10, ones, vstack  
from numpy.linalg import norm, lstsq
from Temp_Schemes.Schemes import *
from ODEs.ODE import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, newton
from Integrador import Integrate_ODE


def Error_ODE( t, order, U0, F, Scheme): # Order = Euler =1, RK4 = 4...
                           
       N = len(t)-1; Nv = len(U0)  # N = número evaluaciones, Nv = Número variables (x,y,vx,vy...)
       t1 = t                      # Tiempo primera malla
       t2 = zeros(2*N+1)           # Tiempo malla refinada
       Error = zeros((Nv, N))
       
       for i in range(N):  
           t2[2*i]   = t1[i] 
           t2[2*i+1] = ( t1[i] + t1[i+1] )/2
       t2[2*N] = t1[N]
      
       
       U1 =   Integrate_ODE(U0, F, t1,  Scheme) 
       U2 =   Integrate_ODE(U0, F, t2,  Scheme)    
       
       for i in range(N):  
            Error[:,i] = ( U2[:, 2*i]- U1[:, i] )/( 1 - 1./2**order ) 
        
       Solution = U1 + Error 
       
       return Error, Solution 
                 
def Temporal_convergence_rate( t, m, U0, F, scheme ): # Numero de mallas
                               
     
      log_E = zeros(m) 
      log_N = zeros(m)
      N = len(t)-1
      t1 =t
      U1 = Integrate_ODE( U0 ,F, t1, scheme) 
     

      for i in range(m): 
         N =  2 * N 
         t2 = array( zeros(N+1) )
         t2[0:N+1:2] =  t1; t2[1:N:2] = ( t1[1:int(N/2)+1]  + t1[0:int(N/2)] )/ 2 
         U2 = Integrate_ODE( U0, F, t2,  scheme)           
        
         error = norm( U2[:, N-1] - U1[:, int( (N-1)/2 )] ) 
         log_E[i] = log10( error );  log_N[i] = log10( N )
         t1 = t2;  U1 = U2;   

     
      for j in range(m): 
         if abs(log_E[j]) > 12 :  break 
      #j = min(j, m-1) 
      x = log_N[0:j+1];  y = log_E[0:j+1]
      A = vstack( [ x, ones(len(x)) ] ).T
      m, c = lstsq(A, y, rcond=None)[0]
      order = abs(m) 
      log_E = log_E - log10( 1 - 1./2**order) 

      return order, log_E, log_N