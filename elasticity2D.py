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
W = 1
mu = 1
rho = 1
delta = W / L
gamma = 0.4 * delta**2
beta = 1.25
lambda_ = beta
g = gamma

# Create mesh and define function space
mesh = RectangleMesh(Point(0, 0), Point(L, W), 4, 4)
mesh = Mesh("geometry.xml")
coords = mesh.coordinates()
x = coords[:, 0]
y = coords[:, 1]
x[:] = (x - min(x)) / (max(x) - min(x)) * W
y[:] = (y - min(y)) / (max(y) - min(y)) * L
V = VectorFunctionSpace(mesh, 'P', 1)

# Mark boundary subdomians
boundary_parts = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
boundary_parts.set_all(0)


class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1, 1e-4)
Gamma_Top = TopBoundary()


class PointLoad(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.5, 1e-4) and near(x[1], 0, 1e-4)
Gamma_Point = PointLoad()
Gamma_Point.mark(boundary_parts, 2)

# bc = DirichletBC(V, Constant((0, 0)), clamped_boundary)
bct = DirichletBC(V, Constant((0.0, 0.0)), Gamma_Top, method="pointwise")
T = Constant((0.0, delta * -2e-1)) # Point load on the boundary
bcp = DirichletBC(V, T, Gamma_Point, method="pointwise")
bcs = [bct, bcp]


# Define strain and stress


def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)
    # return sym(nabla_grad(u))


def sigma(u):
    return lambda_ * nabla_div(u) * Identity(d) + 2 * mu * epsilon(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()  # space dimension
v = TestFunction(V)
# f = Constant((0, -rho * g))
B = Constant((0.0,   0.0)) # Body force per unit volume
dss = ds(subdomain_data=boundary_parts)
a = inner(sigma(u), epsilon(v)) * dx
L = dot(B, v) * dx + dot(T, v) * dss(2)

# Compute solution
A, b = assemble_system(a, L, bcs)
u = Function(V)
U = u.vector()
# solve(a == L, u, b)
solve(A, U, b)

# Plot solution
plot(u, title='Displacement', mode='displacement')

# Plot stress
s = sigma(u) - (1. / 3) * tr(sigma(u)) * Identity(d)  # deviatoric stress
von_Mises = sqrt(3. / 2 * inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)
plot(von_Mises, title='Stress intensity')

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
print('min/max u:',
      u_magnitude.vector().get_local().min(),
      u_magnitude.vector().get_local().max())

# Save solution to file in VTK format
File('elasticity2D/displacement.pvd') << u
File('elasticity2D/von_mises.pvd') << von_Mises
File('elasticity2D/magnitude.pvd') << u_magnitude
