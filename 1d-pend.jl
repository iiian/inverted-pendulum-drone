# A model of the equations of motions for simple 1D inverted pendulum problem

using SymPy

# vars & such
@vars t l g m M
x, y, theta, F = symbols("x y theta F", cls=sympy.Function)

x, y, theta, F = x(t), y(t), theta(t), F(t)

# kinematics
dtx, dty, dttheta = diff(x, t), diff(y, t), diff(theta, t)
x_m = x - l * sin(theta)
y_m = l * cos(theta)

dtx_m, dty_m = diff(x_m, t), diff(y_m, t)

# Energy equations
PE = m * g * y_m

v_M = sqrt(dtx^2 + dty^2)
v_m = sqrt(dtx_m^2 + dty_m^2)
KE = simplify(Rational(1, 2) * (M * v_M^2 + m * v_m^2))

# Lagrangian

L = simplify(KE - PE)

Qi_x = diff(diff(L, dtx), t) - diff(L, x)

# Notice that Qi_y bleeds out to literally just be the good old "F = ma". This is also important, because, we can meaningfully set this Qi_y to 0; There will never be a force applied in the Y direction.
Qi_y = diff(diff(L, dty), t) - diff(L, y)

# Qi_theta ends up being a wild term, but we'll set it equal to zero, as it is a conserved quantity.
Qi_theta = diff(diff(L, dttheta), t) - diff(L, theta)