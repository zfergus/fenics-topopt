from __future__ import print_function, division
from fenics import *
import mshr

from elasticity.boundary_conditions import BoundaryConditions
from elasticity.utils import scale_mesh
from run_simulation import run_simulation


class LBracketBoundaryConditions(BoundaryConditions):
    def get_fixed(self):
        width, height, tol = self.width, self.height, self.tol

        class TopBoundary(SubDomain):
            """ Constrain the bottom to not move. """
            def inside(self, x, on_boundary):
                return near(x[1], height, tol)
        return [TopBoundary()]

    def get_forces(self):
        width, height, tol = self.width, self.height, self.tol

        class PointLoad(SubDomain):
            """ Add a point load to the top center. """
            def inside(self, x, on_boundary):
                return near(x[0], width, tol) and near(x[1], height / 3., tol)
        return [PointLoad()], [Constant((0.0, -3e-1))]


if __name__ == "__main__":
    bc = LBracketBoundaryConditions(1, 1, 5e-2)

    # Create the mesh to solve linear elasticity on.
    large_quad = mshr.Polygon([Point(0, 0), Point(1, 0), Point(1, 1),
        Point(0, 1)])
    small_quad = mshr.Polygon([Point(1 / 3., 1 / 3.), Point(1, 1 / 3.),
        Point(1, 1), Point(1 / 3., 1)])
    domain = large_quad - small_quad
    mesh = mshr.generate_mesh(domain, 30)
    run_simulation(mesh, bc, "L-bracket/v100-")

    mesh = Mesh("meshes/L-bracket-v20.xml")
    scale_mesh(mesh, 1, 1)
    run_simulation(mesh, bc, "L-bracket/v20-")

    mesh = Mesh("meshes/L-bracket-v10.xml")
    scale_mesh(mesh, 1, 1)
    run_simulation(mesh, bc, "L-bracket/v10-")
