from __future__ import print_function, division
from fenics import *
import mshr

from elasticity.utils import scale_mesh
from run_simulation import run_simulation


class BoundaryConditions:
    def __init__(self, width, height, tol):
        self.width, self.height, self.tol = width, height, tol

    def get_fixed(self):
        width, height, tol = self.width, self.height, self.tol

        class BottomSupport(SubDomain):
            """ Add a point load to the top center. """
            def inside(self, x, on_boundary):
                return near(x[1], -1, tol)
        return [BottomSupport()]

    def get_forces(self):
        width, height, tol = self.width, self.height, self.tol

        class TopLoad(SubDomain):
            """ Constrain the bottom to not move. """
            def inside(self, x, on_boundary):
                return near(x[0], 0, 0.25) and near(x[1], 1, tol)
        return [TopLoad()], [Constant((0.0, -1))]

    def get_body_force(self):
        return Constant((0, 0))


if __name__ == "__main__":
    width, height, tol, E = 2, 2, 5e-2, 5e1
    mesh = Mesh("meshes/circle-in-square.xml")
    bc = BoundaryConditions(width, height, tol)
    run_simulation(mesh, bc, "circle/", E=E)
