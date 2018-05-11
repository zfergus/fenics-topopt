from __future__ import print_function, division
from fenics import *
import mshr

from linear_elasticity import linear_elasticity
from von_Mises_stress import von_Mises_stress
from utils import scale_mesh


def main():
    height = 1
    width = 1
    delta = width / height

    # Create the mesh to solve linear elasticity on.
    large_quad = mshr.Polygon([Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)])
    small_quad = mshr.Polygon([Point(0.5, 0.5), Point(1, 0.5), Point(1, 1), Point(0.5, 1)])
    domain = large_quad - small_quad
    mesh = mshr.generate_mesh(domain, 30)

    scale_mesh(mesh, width, height)

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
            return near(x[0], width, tol) and near(x[1], height / 2., tol)
    gamma_point = PointLoad()
    gamma_point.mark(boundary_parts, 2)

    B = Constant((0.0, 0.0)) # Body force per unit volume
    bct = DirichletBC(V, Constant((0.0, 0.0)), gamma_top, method="pointwise")

    nframes = 30 * 10
    file = File("output/L-bracket/L-bracket.pvd")
    for t in range(nframes):
        T = Constant((0.0, -(t + 1) / float(nframes))) # Point load on the boundary

        # Boundary conditions on the subdomains
        bcp = DirichletBC(V, T, gamma_point, method="pointwise")
        bcs = [bct, bcp]

        dss = ds(subdomain_data=boundary_parts)
        L = lambda v: dot(B, v) * dx + dot(T, v) * dss(2)

        u = linear_elasticity(V, L, bcs)
        file << u
        print("{:d}/{:d}".format(t, nframes))

    # Compute magnitude of displacement
    # V = FunctionSpace(mesh, "P", 1)
    # u_magnitude = sqrt(dot(u, u))
    # u_magnitude = project(u_magnitude, V)
    # file << u_magnitude
    #
    # print("min/max u: {:g}, {:g}".format(
    #     u_magnitude.vector().get_local().min(),
    #     u_magnitude.vector().get_local().max()))
    #
    # # Save solution to file in VTK format
    # file << von_Mises_stress(mesh, u)

if __name__ == "__main__":
    main()
