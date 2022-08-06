# A model of the equations of motions for simple 1D inverted pendulum problem

using SymPy

# vars & such
@vars t l_alpha l_beta m1 m2 g
alpha, beta = symbols("alpha beta", cls=sympy.Function)
alpha, beta = alpha(t), beta(t)

dalpha, dbeta = alpha.diff(t), beta.diff(t)

# kinematics
x1 = l_alpha * sin(alpha)
y1 = -l_alpha * cos(alpha)

x2 = x1 + l_beta * sin(beta)
y2 = y1 - l_beta * cos(beta)

# naively derived velocities from kinematics
# -- am I missing something?
dx1, dx2, dy1, dy2 = x1.diff(t), x2.diff(t), y1.diff(t), y2.diff(t)

v1, v2 = sqrt(dx1^2 + dy1^2), sqrt(dx2^2 + dy2^2)

# Energy equations

PE = m1 * g * y1 + m2 * g * y2
KE = Rational(1, 2) * m1 * v1^2 + Rational(1, 2) * m2 * v2^2

# Lagrangian
L = KE - PE

dLT(L, t, dqi) = L.diff(dqi, t)
dLS(L, qi) = L.diff(qi)

# wrt alpha
dLT_alpha = L.diff(dalpha, t)
dLS_alpha = L.diff(alpha)

dLT_beta = L.diff(dbeta, t)
dLS_beta = L.diff(beta)

# solve for state

print(sympy.latex(Eq(dLT_alpha - dLS_alpha, 0)))