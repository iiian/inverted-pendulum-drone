# Derived from https://www.flyingmachinearena.ethz.ch/wp-content/publications/2011/hehn_dandrea_flying_inverted_pendulum.pdf
using LinearAlgebra, ControlSystems

#=
  Compute FSFB Controller K matrix based on
  1. linearized system dynamics derived from taylor series approximation. 
  2. dx/dt = (A - BK)x, where K is the LQR solution given A, B, and Q/R matrices

  x/y & z systems are separated

  below are three lqr computations for three separated flight controllers based on the total 13 state, 4 control system for the drone:
    an x-axis controller
    a y-axis controller
    and a z-axis controller

  The x & y axis controllers are both 5-state single input controllers: 
    states:
      x/y, 
      resp. velocity, 
      resp. euler angle of vehicle body, 
      resp. pendulum-COM offset diff component, 
      its resp. velocity
    control:
      x/y rotational body rate (ω x/y)

  The z axis controller is a 2-state single input controller:
    states:
      z position of vehicle body
      z component velocity of vehicle body

    control:
      acceleration produced by quadcopter thrust
=#


# Tunable cost paramaters for LQR. 
#   Q matrix is set up to penalize cost of position of drone
#   R matrix is set up to penalize cost of control inputs (ie battery life, want to move optimally, but quickly)
qQ = [
  1 0 0 0 0
  0 0 0 0 0
  0 0 0 0 0
  0 0 0 0 0
  0 0 0 0 0
]

qR = 100 * I

# Produce the A matrix for the x & y controllers based on parameters
# g : gravitational constant
# L : length of pendulum to center of mass
# isX: true if for x controller, false for y controller
function qA(g, L, isX::Bool)
  sign = isX ? 1 : -1
  [
    0 1 0 0 0
    0 0 sign*g 0 0
    0 0 0 0 0
    0 0 0 0 1
    0 0 -sign*g g/L 0
  ]
end

g, L = 9.81, 1
xA = qA(g, L, true)
yA = qA(g, L, false)

# identical for both x and y controllers 
qB = [
  0
  0
  1
  0
  0
]

zQ = [
  1 0
  0 0
]

zR = qR

zA = [
  0 1
  0 0
]

zB = [0; 1]

xK = lqr(Continuous, xA, qB, qQ, qR)
yK = lqr(Continuous, yA, qB, qQ, qR)
zK = lqr(Continuous, zA, zB, zQ, zR)

# Simulation:
#  Given some simple initial condition, we can update the system continuously and observe it approaching the zero which the system has been linearized about.

x = [
  10 # x coord
  0.5 # x speed
  20 * 2 * pi / 360 # euler angle (rads) 
  0.1 # center of mass offset
  0.2 # COM offset speed 
]

y = [
  4 # y coord
  2 # y speed
  20 * 2 * pi / 360 # euler angle (rads)
  0.1 # center of mass offset (y)
  0.06 # COM offset speed
]

z = [
  5 # z coord
  0 # z speed
]

print("Initial condition for x state is: ")
println(x)
print("Initial condition for y state is: ")
println(y)
print("Initial condition for z state is: ")
println(z)

# iterate for one thousand 100ms time steps steps
for i in 1:1000
  global x += (xA - (qB * xK)) * x * 0.1
  global y += (yA - (qB * yK)) * y * 0.1
  global z += (zA - (zB * zK)) * z * 0.1
end

print("Final condition for x state is: ")
println(x)
print("Final condition for y state is: ")
println(y)
print("Final condition for z state is: ")
println(z)