from __future__ import print_function, division
from fenics import *
import mshr

from elasticity.linear_elasticity import linear_elasticity
from elasticity.von_Mises_stress import von_Mises_stress
from elasticity.utils import scale_mesh


def main():
    height = 100
    width = 600
    delta = width / height

    # Create the mesh to solve linear elasticity on.
    mesh = Mesh("meshes/bridge2.xml")

    scale_mesh(mesh, width, height)

    # Define the function space
    V = VectorFunctionSpace(mesh, "P", 1)
    # Mark boundary subdomians
    boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundary_parts.set_all(0)

    # Tolarance of boundary near checks.
    tol = 5e-2

    class BottomBoundary(SubDomain):
        """ Constrain the bottom to not move. """
        def inside(self, x, on_boundary):
            return ((near(x[0], 0.1, 2e1) or near(x[0], width - 0.1, 2e1)) and
                near(x[1], 0, tol))
    gamma_bottom = BottomBoundary()

    class PointLoad(SubDomain):
        """ Add a point load to the top center. """
        def inside(self, x, on_boundary):
            return (near(x[0], width / 2.0, width / 4) and
                near(x[1], height, tol))
    gamma_point = PointLoad()
    gamma_point.mark(boundary_parts, 2)

    B = Constant((0.0, 0.0)) # Body force per unit volume
    T = Constant((0.0, -2e0)) # Point load on the boundary

    # Boundary conditions on the subdomains
    bct = DirichletBC(V, Constant((0.0, 0.0)), gamma_bottom)
    bcs = [bct]

    dss = ds(domain=mesh, subdomain_data=boundary_parts)
    L = lambda v: dot(B, v) * dx + dot(T, v) * dss(2)

    u = linear_elasticity(V, L, bcs, E=1e3)
    File("output/MBB/displacement.pvd") << u

    # Compute magnitude of displacement
    V = FunctionSpace(mesh, "P", 1)
    u_magnitude = sqrt(dot(u, u))
    u_magnitude = project(u_magnitude, V)
    File("output/MBB/magnitude.pvd") << u_magnitude

    print("min/max u: {:g}, {:g}".format(
        u_magnitude.vector().get_local().min(),
        u_magnitude.vector().get_local().max()))

    # Save solution to file in VTK format
    File("output/MBB/von_mises.pvd") << von_Mises_stress(mesh, u)

if __name__ == "__main__":
    main()
