from __future__ import print_function, division
from fenics import *
import collections


def epsilon(u):
    """ Define strain. """
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)


def sigma(u, lmda=1.25, mu=1):
    """ Define stress. """
    d = u.geometric_dimension()
    return lmda * nabla_div(u) * Identity(d) + 2 * mu * epsilon(u)


def linear_elasticity(V, B, Ts, bcs, boundary_parts):
    """
    Solve the Linear elatic problemu using FEniCS:
        -div(sigma(u)) = f
    Solves the problem on the given mesh with the given MeshFunction for the
    boundary conditions.
    """
    # Define variational problem
    u = TrialFunction(V)
    d = u.geometric_dimension() # space dimension
    v = TestFunction(V)
    dss = ds(subdomain_data=boundary_parts)
    a = inner(sigma(u), epsilon(v)) * dx
    if isinstance(Ts, Coefficient):
        Ts = [Ts]
    T = sum([dot(T, v) * dss(i) for i, T in enumerate(Ts)])
    L = dot(B, v) * dx + T
    # Compute solution
    A, b = assemble_system(a, L, bcs)
    u = Function(V)
    U = u.vector()

    # solve(a == L, u, b)
    solve(A, U, b)
    return u
