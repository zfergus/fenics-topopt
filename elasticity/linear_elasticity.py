from __future__ import print_function, division
from fenics import *
import collections


def epsilon(u):
    """ Define strain. """
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)


def sigma(u, mu=1, lmbda=1.25):
    """ Define stress. """
    d = u.geometric_dimension()
    return lmbda * nabla_div(u) * Identity(d) + 2 * mu * epsilon(u)


def linear_elasticity(V, L, bcs, E=1, nu=0.3):
    """
    Solve the Linear elatic problemu using FEniCS:
        -div(sigma(u)) = f
    Solves the problem on the given mesh with the given MeshFunction for the
    boundary conditions.
    """
    mu, lmbda = E / (2. * (1 + nu)), E * nu / ((1. + nu) * (1. - 2. * nu))
    # Define variational problem
    u = TrialFunction(V)
    d = u.geometric_dimension() # space dimension
    v = TestFunction(V)
    a = inner(sigma(u, mu, lmbda), epsilon(v)) * dx
    # Compute solution
    L = L(v)
    u = Function(V)
    solve(a == L, u, bcs)
    return u
