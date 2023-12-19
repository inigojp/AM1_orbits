# AM1_orbits
Iñigo Javier Palacios Martínez
El GitHub está estructurado de la siguiente forma:

- Esquemas temporales en Carpeta Temp_Schemes, archivo Schemes.py, RK4, Crank-Nicolson, Euler y Euler inverso.
- Ecuaciones diferenciales en ODES, ODE.py. Órbitas keplerianas, problema de los N cuerpos, Osciladores...
- Funciones para determinar el error numérico mediante extrapolación de Richardson y para calcular la convergencia de los métodos. En carpeta Error, Richardson.py y Stability.py

Dos funciones para integrar las ecuaciones en Integrador.py

- Hito 1: Integrar órbita de Kepler con métodos RK4, Crank-Nicolson y Euler
- Hito 2: Separar métodos en modulos, integrar problema genérico de Cauchy
- Hito 3: Estimación de error mediante interpolacion de Richardson
- Hito 4: Regiones de estabilidad 
- Hito 5: Integración del problema de los N cuerpos
- Hito 6: Integración CR3BP, puntos de Lagrange
- Hito 7: Integrar órbita de Arenstof mediante GBS, RK4 y RK4 de paso variable