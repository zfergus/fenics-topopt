from __future__ import print_function, division
from fenics import *
import mshr

from linear_elasticity import linear_elasticity
from von_Mises_stress import von_Mises_stress
from utils import scale_mesh


def main():
    # Create the mesh to solve linear elasticity on.
    mesh = Mesh("meshes/circle-in-square.xml")

    # Define the function space
    V = VectorFunctionSpace(mesh, "P", 1)
    # Mark boundary subdomians
    boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary_parts.set_all(0)

    # Tolarance of boundary near checks.
    tol = 2e-4

    class TopLoad(SubDomain):
        """ Constrain the bottom to not move. """
        def inside(self, x, on_boundary):
            return near(x[0], 0, 0.2) and near(x[1], 1, 1e-3)
    gamma_top = TopLoad()
    gamma_top.mark(boundary_parts, 1)

    class BottomLoad(SubDomain):
        """ Add a point load to the top center. """
        def inside(self, x, on_boundary):
            return near(x[0], 0, 0.2) and near(x[1], -1, 1e-3)
    gamma_bottom = BottomLoad()
    gamma_bottom.mark(boundary_parts, 2)

    B = Constant((0.0, 0.0)) # Body force per unit volume
    T1 = Constant((0.0, -1e-2)) # Point load on the top boundary
    T2 = Constant((0.0, 1e-2)) # Point load on the bottom boundary

    # Boundary conditions on the subdomains
    bct = DirichletBC(V, T1, gamma_top, method="topological")
    bcb = DirichletBC(V, T2, gamma_bottom, method="topological")
    bcs = [bct, bcb]

    u = linear_elasticity(V, B, [T1, T2], bcs, boundary_parts)
    File("output/circle/displacement.pvd") << u

    # Compute magnitude of displacement
    V = FunctionSpace(mesh, "P", 1)
    u_magnitude = sqrt(dot(u, u))
    u_magnitude = project(u_magnitude, V)
    File("output/circle/magnitude.pvd") << u_magnitude

    print("min/max u: {:g}, {:g}".format(
        u_magnitude.vector().get_local().min(),
        u_magnitude.vector().get_local().max()))

    # Save solution to file in VTK format
    File("output/circle/von_mises.pvd") << von_Mises_stress(mesh, u)

if __name__ == "__main__":
    main()
