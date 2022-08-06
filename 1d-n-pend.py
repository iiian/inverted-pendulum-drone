# based off of https://www.moorepants.info/blog/npendulum.html

from sympy import symbols
import sympy
from sympy.physics.mechanics import *

n = 5

q = dynamicsymbols('q:' + str(n + 1)) #generalized coordinates
u = dynamicsymbols('u:' + str(n + 1)) #generalized speeds
f = dynamicsymbols('f') # force applied

m = symbols('m:' + str(n + 1)) # mass of each bob
l = symbols('l:' + str(n)) # length of ea. link
g, t = symbols('g t')

I = ReferenceFrame('I') # inertial reference frame
O = Point('O') # origin
O.set_vel(I, 0) # origin velocity is zero

# Define the first point of the pendulum as a particle which has mass, and can move laterally:
P0 = Point('P0')
P0.set_pos(O, q[0] * I.x) # set the position of P0
# ^ this says: wrt to origin, P0 has point q_0 along the x axis
P0.set_vel(I, u[0] * I.x)
Pa0 = Particle('Pa0', P0, m[0])

frames = [I]
points = [P0]
particles = [Pa0]
forces = [(P0, f * I.x - m[0] * g * I.y)]
kindiffs = [q[0].diff(t) - u[0]]

for i in range(n):
    Bi = I.orientnew('B' + str(i), 'Axis', [q[i + 1], I.z])   # Create a new frame
    Bi.set_ang_vel(I, u[i + 1] * I.z)                         # Set angular velocity
    frames.append(Bi)                                         # Add it to the frames list

    Pi = points[-1].locatenew('P' + str(i + 1), l[i] * Bi.x)  # Create a new point
    Pi.v2pt_theory(points[-1], I, Bi)                         # Set the velocity
    points.append(Pi)                                         # Add it to the points list

    Pai = Particle('Pa' + str(i + 1), Pi, m[i + 1])           # Create a new particle
    particles.append(Pai)                                     # Add it to the particles list

    forces.append((Pi, -m[i + 1] * g * I.y))                  # Set the force applied at the point

    kindiffs.append(q[i + 1].diff(t) - u[i + 1])              # Define the kinematic ODE:  dq_i / dt - u_i = 0

kane = KanesMethod(I, q_ind=q, u_ind=u, kd_eqs=kindiffs) # Initialize the object
fr = kane.kanes_equations(forces, particles)     # Generate EoM's fr + frstar = 0

print(sympy.latex(fr))
