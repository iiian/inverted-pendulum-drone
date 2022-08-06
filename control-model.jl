using SymPy
# time
@vars t
# define parameteric state variables
# 10D state
x1, y1, z1, theta, psi = symbols("x1 y1 z1 theta psi", cls=sympy.Function)
x1, y1, z1, theta, psi = x1(t), y1(t), z1(t), theta(t), psi(t)
dx1, dy1, dz1, dtheta, dpsi = diff(x1, t), diff(y1, t), diff(z1, t), diff(theta, t), diff(psi, t)
# static parts of eqn
@vars M m g r

# build eqn constraints
x2 = r * (dtheta * cos(theta) * cos(psi) - dpsi * sin(theta) * sin(psi)) + dx1
y2 = r * (dtheta * cos(theta) * sin(psi) + dpsi * sin(theta) * cos(psi)) + dy1
z2 = r * diff(theta, t) * sin(theta) + z1
dx2 = diff(x2, t)
dy2 = diff(y2, t)
dz2 = diff(z2, t)

# build lagrangian
PE = (M * g * z1) + (m * g * (z1 + r * cos(theta)))
v1 = sqrt(dx1^2 + dy1^2 + dz1^2)
v2 = sqrt(dx2^2 + dy2^2 + dz2^2)
KE = 0.5 * M * v1^2 + 0.5 * m * v2^2
L = KE - PE

# TODO:

# compute jacobian

# get 10D linear, find optimal eigenstate reshape using LQR
