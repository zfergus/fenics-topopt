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


def linear_elasticity(V, L, bcs):
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
    a = inner(sigma(u), epsilon(v)) * dx
    # Compute solution
    L = L(v)
    A, b = assemble_system(a, L, bcs)
    u = Function(V)
    U = u.vector()

    # solve(a == L, u, b)
    solve(A, U, b)
    return u
