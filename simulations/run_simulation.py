from fenics import *

from elasticity.linear_elasticity import linear_elasticity
from elasticity.von_Mises_stress import von_Mises_stress


def run_simulation(mesh, bc, filename_prefix, E=1e2, nu=0.3):
    u = linear_elasticity(mesh, bc, E=E)

    # Compute magnitude of displacement
    V = FunctionSpace(mesh, "P", 1)
    u_magnitude = sqrt(dot(u, u))
    u_magnitude = project(u_magnitude, V)

    print("min/max u: {:g}, {:g}".format(
        u_magnitude.vector().get_local().min(),
        u_magnitude.vector().get_local().max()))

    von_Mises = von_Mises_stress(mesh, u)

    # Save solution to file in VTK format
    File("output/{:s}displacement.pvd".format(filename_prefix)) << u
    File("output/{:s}magnitude.pvd".format(filename_prefix)) << u_magnitude
    File("output/{:s}von_mises.pvd".format(filename_prefix)) << von_Mises
