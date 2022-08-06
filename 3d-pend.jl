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
x2 = r * (dtheta * cos(theta) * cos(psi) - dpsi * sin(theta) * sin(psi)) + dx1 #TODO: could this all be easier in spherical coordinates?
y2 = r * (dtheta * cos(theta) * sin(psi) + dpsi * sin(theta) * cos(psi)) + dy1
z2 = r * sin(psi) + z1
dx2 = diff(x2, t)
dy2 = diff(y2, t)
dz2 = diff(z2, t)

# build lagrangian
PE = (M * g * z1) + (m * g * (z1 + r * cos(theta)))
v1 = sqrt(dx1^2 + dy1^2 + dz1^2)
v2 = sqrt(dx2^2 + dy2^2 + dz2^2)
KE = Rational(1, 2) * M * v1^2 + Rational(1, 2) * m * v2^2
L = KE - PE

# L differentiated wrt acceleration of system
dL_dtheta_dt = diff(diff(L, dtheta), t)
dL_dpsi_dt = diff(diff(L, dpsi), t)
dL_dx1_dt = diff(diff(L, dx1), t)
dL_dy1_dt = diff(diff(L, dy1), t)
dL_dz1_dt = diff(diff(L, dz1), t)

# L differentiated wrt to state of system
dL_theta = diff(L, theta)
dL_psi = diff(L, psi)
dL_x1 = diff(L, x1)
dL_y1 = diff(L, y1)
dL_z1 = diff(L, z1)


Qi_theta = simplify(dL_dtheta_dt - dL_theta)
Qi_psi = simplify(dL_dpsi_dt - dL_psi)
Qi_x1 = simplify(dL_dx1_dt - dL_x1)
Qi_y1 = simplify(dL_dy1_dt - dL_y1)
Qi_z1 = simplify(dL_dz1_dt - dL_z1)

print("here goes...")
theta_accel = solve(Eq(Qi_theta, 0), theta.diff(t, t))
print("made it")
# At this point, we could try and compute out the 10D state vector's rate of change:
# psi_accel = solve(Eq(Qi_psi, 0), psi.diff(t, t))
# x1_accel = solve(Eq(Qi_x1, 0), x1.diff(t, t))
# y1_accel = solve(Eq(Qi_y1, 0), y1.diff(t, t))
# z1_accel = solve(Eq(Qi_z1, 0), z1.diff(t, t))

# theta_vel = solve(Eq(Qi_theta, 0), theta.diff(t))
# psi_vel = solve(Eq(Qi_psi, 0), psi.diff(t))
# x1_vel = solve(Eq(Qi_x1, 0), x1.diff(t))
# y1_vel = solve(Eq(Qi_y1, 0), y1.diff(t))
# z1_vel = solve(Eq(Qi_z1, 0), z1.diff(t))

# The problem is that the equations are prohibitively absurd and complex.
# After reading through https://www.flyingmachinearena.ethz.ch/wp-content/publications/2011/hehn_dandrea_flying_inverted_pendulum.pdf
# it seems like the problem comes from defining the attitude of the pendulum in terms of just two angles: theta and psi.
# If, instead, I were to use euler angles to first create a transformation between inertial frames, this shakes out to yield a state vector in 7 variables, which although complex
# seems much less complex than my naive attempt, and is smaller than my 10D attempt.