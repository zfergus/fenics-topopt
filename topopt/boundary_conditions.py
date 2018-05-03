from __future__ import division

import numpy as np


class BoundaryConditions(object):

    def __init__(self, nelx, nely):
        self.nelx, self.nely = nelx, nely

    def get_fixed_nodes(self):
        # Return a list of fixed nodes for the problem
        dofs = np.arange(2 * (self.nelx + 1) * (self.nely + 1))
        fixed = np.union1d(dofs[0:2 * (self.nely + 1):2],
            np.array([2 * (self.nelx + 1) * (self.nely + 1) - 1]))
        return fixed

    def get_forces(self):
        # Return the force vector for the problem
        ndof = 2 * (self.nelx + 1) * (self.nely + 1)
        f = np.zeros((ndof, 1))
        f[1, 0] = -1
        return f

    def get_passive_elements(self):
        return
