from __future__ import print_function, division
from fenics import *
import mshr

from elasticity.linear_elasticity import linear_elasticity
from elasticity.von_Mises_stress import von_Mises_stress
from elasticity.utils import scale_mesh


def main(mesh, filename_prefix):
    width, height = 1, 1
    # Define the function space
    V = VectorFunctionSpace(mesh, "P", 1)
    # Mark boundary subdomians
    boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary_parts.set_all(0)

    # Tolarance of boundary near checks.
    tol = 2e-4

    class TopBoundary(SubDomain):
        """ Constrain the bottom to not move. """
        def inside(self, x, on_boundary):
            return near(x[1], height, tol)
    gamma_top = TopBoundary()

    class PointLoad(SubDomain):
        """ Add a point load to the top center. """
        def inside(self, x, on_boundary):
            return near(x[0], width, tol) and near(x[1], height / 3., tol)
    gamma_point = PointLoad()
    gamma_point.mark(boundary_parts, 2)

    B = Constant((0.0, 0.0)) # Body force per unit volume
    T = Constant((0.0, -2e-1)) # Point load on the boundary

    # Boundary conditions on the subdomains
    bct = DirichletBC(V, Constant((0.0, 0.0)), gamma_top, method="pointwise")
    bcp = DirichletBC(V, T, gamma_point, method="pointwise")
    bcs = [bct, bcp]

    dss = ds(subdomain_data=boundary_parts)
    L = lambda v: dot(B, v) * dx + dot(T, v) * dss(2)

    u = linear_elasticity(V, L, bcs)

    # Compute magnitude of displacement
    V = FunctionSpace(mesh, "P", 1)
    u_magnitude = sqrt(dot(u, u))
    u_magnitude = project(u_magnitude, V)

    print("min/max u: {:g}, {:g}".format(
        u_magnitude.vector().get_local().min(),
        u_magnitude.vector().get_local().max()))

    von_Mises = von_Mises_stress(mesh, u)

    # Save solution to file in VTK format
    File("output/L-bracket/{:s}displacement.pvd".format(filename_prefix)) << u
    (File("output/L-bracket/{:s}magnitude.pvd".format(filename_prefix)) <<
        u_magnitude)
    (File("output/L-bracket/{:s}von_mises.pvd".format(filename_prefix)) <<
        von_Mises)

if __name__ == "__main__":
    # Create the mesh to solve linear elasticity on.
    large_quad = mshr.Polygon([Point(0, 0), Point(1, 0), Point(1, 1),
        Point(0, 1)])
    small_quad = mshr.Polygon([Point(1 / 3., 1 / 3.), Point(1, 1 / 3.),
        Point(1, 1), Point(1 / 3., 1)])
    domain = large_quad - small_quad
    mesh = mshr.generate_mesh(domain, 30)
    main(mesh, "")
