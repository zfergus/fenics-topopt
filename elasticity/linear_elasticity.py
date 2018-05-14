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


def linear_elasticity(mesh, bc, E=1, nu=0.3):
    """
    Solve the Linear elatic problemu using FEniCS:
        -div(sigma(u)) = f
    Solves the problem on the given mesh with the given MeshFunction for the
    boundary conditions.
    """
    # Define the function space
    V = VectorFunctionSpace(mesh, "P", 1)
    # Mark boundary subdomians
    boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary_parts.set_all(0)

    # Get fixed subdomains
    fixed_subdomains = bc.get_fixed()
    # Enforce the fixed subdomain as a Dirichlet boundary with value 0.
    fixed_bcs = [DirichletBC(V, Constant((0.0, 0.0)), fixed)
        for fixed in fixed_subdomains]

    # Mark the force dubdomains
    forces, Ts = bc.get_forces()
    for i, force in enumerate(forces):
        force.mark(boundary_parts, len(fixed_subdomains) + i)
    B = bc.get_body_force() # Body force per unit volume

    # Define variational problem
    u = TrialFunction(V)
    d = u.geometric_dimension() # space dimension
    v = TestFunction(V)
    mu, lmbda = E / (2. * (1 + nu)), E * nu / ((1. + nu) * (1. - 2. * nu))
    a = inner(sigma(u, mu, lmbda), epsilon(v)) * dx
    # Boundary conditions on the subdomains
    dss = ds(domain=mesh, subdomain_data=boundary_parts)
    L = dot(B, v) * dx + sum([dot(T, v) * dss(len(fixed_subdomains) + i)
        for i, T in enumerate(Ts)])
    # Compute solution
    u = Function(V)
    solve(a == L, u, fixed_bcs)
    return u
