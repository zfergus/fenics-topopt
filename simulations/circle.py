from __future__ import print_function, division
from fenics import *
import mshr

from elasticity.linear_elasticity import linear_elasticity
from elasticity.von_Mises_stress import von_Mises_stress
from elasticity.utils import scale_mesh


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
            return near(x[0], 0, 0.25) and near(x[1], 1, 1e-3)
    gamma_top = TopLoad()
    gamma_top.mark(boundary_parts, 1)

    class BottomSupport(SubDomain):
        """ Add a point load to the top center. """
        def inside(self, x, on_boundary):
            return near(x[1], -1, 1e-3)
    gamma_bottom = BottomSupport()

    B = Constant((0.0, 0.0)) # Body force per unit volume
    T = Constant((0.0, 0)) # Point load on the top boundary

    # Boundary conditions on the subdomains
    bct = DirichletBC(V, T, gamma_top, method="pointwise")
    bcb = DirichletBC(V, Constant((0, 0)), gamma_bottom, method="pointwise")
    bcs = [bct, bcb]

    dss = ds(subdomain_data=boundary_parts)
    L = lambda v: dot(B, v) * dx + dot(T, v) * dss(2)
    u = linear_elasticity(V, L, bcs)
    File("output/circle/displacement-before.pvd") << u

    T = Constant((0.0, -1e-1)) # Point load on the top boundary
    # Boundary conditions on the subdomains
    bct = DirichletBC(V, T, gamma_top, method="pointwise")
    bcb = DirichletBC(V, Constant((0, 0)), gamma_bottom, method="pointwise")
    bcs = [bct, bcb]

    dss = ds(subdomain_data=boundary_parts)
    L = lambda v: dot(B, v) * dx + dot(T, v) * dss(2)
    u = linear_elasticity(V, L, bcs)
    File("output/circle/displacement-after.pvd") << u

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
