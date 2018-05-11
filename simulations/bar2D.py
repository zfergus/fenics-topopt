"""
FEniCS tutorial demo program: Linear elastic problem.

  -div(sigma(u)) = f

The model is used to simulate an elastic beam clamped at
its left end and deformed under its own weight.
"""

from __future__ import print_function
from fenics import *
import matplotlib.pyplot as plt

# Scaled variables
L = 1
W = 0.2
mu = 1
rho = 1
delta = W / L
gamma = 0.4 * delta**2
beta = 1.25
lambda_ = beta
g = gamma

# Create mesh and define function space
mesh = RectangleMesh(Point(0, 0), Point(L, W), 10, 3)
V = VectorFunctionSpace(mesh, 'P', 1)
tol = 1e-2

def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < tol
bc = DirichletBC(V, Constant((0, 0)), clamped_boundary)

# Define strain and stress
epsilon = lambda u: 0.5 * (nabla_grad(u) + nabla_grad(u).T)
sigma = lambda u: lambda_ * nabla_div(u) * Identity(d) + 2 * mu * epsilon(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()  # space dimension
v = TestFunction(V)
f = Constant((0, -rho * g))
T = Constant((0, 0))
a = inner(sigma(u), epsilon(v)) * dx
L = dot(f, v) * dx + dot(T, v) * ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)

plot(u, title='Displacement', mode='displacement')

# Plot stress
s = sigma(u) - (1. / 3) * tr(sigma(u)) * Identity(d)  # deviatoric stress
von_Mises = sqrt(3. / 2 * inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
print('min/max u:',
      u_magnitude.vector().get_local().min(),
      u_magnitude.vector().get_local().max())

# Save solution to file in VTK format
File('output/bar2D/displacement.pvd') << u
File('output/bar2D/von_mises.pvd') << von_Mises
File('output/bar2D/magnitude.pvd') << u_magnitude

plt.show()
