from __future__ import division
from fenics import *


class BoundaryConditions:
    def __init__(self, width, height, tol):
        self.width, self.height, self.tol = width, height, tol

    def get_fixed(self):
        class EmptySubDomain(SubDomain):
            def inside(self, x, on_boundary):
                return False
        return [EmptySubDomain()]

    def get_forces(self):
        class NullLoad(SubDomain):
            def inside(self, x, on_boundary):
                return False
        return [NullLoad()], [Constant((0, 0))]

    def get_body_force(self):
        return Constant((0, 0))
