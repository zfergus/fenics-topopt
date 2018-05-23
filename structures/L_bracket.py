#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import pdb

import numpy as np

from topopt.gui import GUI
from topopt.boundary_conditions import BoundaryConditions
from topopt.solver import Solver
from topopt.utils import xy_to_id, id_to_xy

# Python 2/3 compatibility
try:
    input = raw_input
except NameError:
    pass


class LessTolerantSolver(Solver):
    def __init__(self, nelx, nely, volfrac, penal, rmin, ft, gui, bc):
        super(LessTolerantSolver, self).__init__(nelx, nely, volfrac, penal,
            rmin, ft, gui, bc)
        # set stopping criteria
        self.opt.set_maxeval(4000)
        self.opt.set_ftol_rel(0.0001)


class LBracketBoundaryConditions(BoundaryConditions):
    def __init__(self, nelx, nely, minx, maxy):
        super(LBracketBoundaryConditions, self).__init__(nelx, nely)
        (self.passive_min_x, self.passive_min_y, self.passive_max_x,
            self.passive_max_y) = [minx, 0, nelx - 1, maxy]

    def get_fixed_nodes(self):
        """ Return a list of fixed nodes for the problem. """
        x = np.arange(self.passive_min_x)
        topx_to_id = np.vectorize(
            lambda x: xy_to_id(x, 0, self.nelx, self.nely))
        ids = topx_to_id(x)
        fixed = np.union1d(2 * ids, 2 * ids + 1)
        return fixed

    def get_forces(self):
        """ Return the force vector for the problem. """
        ndof = 2 * (self.nelx + 1) * (self.nely + 1)
        f = np.zeros((ndof, 1))
        fx = self.nelx
        # fy = (self.nely - self.passive_max_y) // 2 + self.passive_max_y
        for i in range(1, 2):
            fy = self.passive_max_y - 1 + 2 * i
            id = xy_to_id(fx, fy, self.nelx, self.nely)
            f[2 * id + 1, 0] = -1
        return f

    def get_passive_elements(self):
        X, Y = np.mgrid[self.passive_min_x:self.passive_max_x + 1,
            self.passive_min_y:self.passive_max_y]
        pairs = np.vstack([X.ravel(), Y.ravel()]).T
        passive_to_ids = np.vectorize(lambda pair: xy_to_id(*pair,
            nelx=self.nelx - 1, nely=self.nely - 1), signature="(m)->()")
        return passive_to_ids(pairs)


def main(nelx, nely, volfrac, penal, rmin, ft):
    print("Minimum compliance problem with MMA")
    print("ndes: {:d} x {:d}".format(nelx, nely))
    print("volfrac: {:g}, rmin: {:g}, penal: {:g}".format(volfrac, rmin, penal))
    print("Filter method: " + ["Sensitivity based", "Density based"][ft])

    # Allocate design variables (as array), initialize and allocate sens.
    x = volfrac * np.ones(nely * nelx, dtype=float)

    title = "ndes: {:d} x {:d}\nvolfrac: {:g}, rmin: {:g}, penal: {:g}".format(
        nelx, nely, volfrac, rmin, penal)
    gui    = GUI(nelx, nely, title)
    solver = LessTolerantSolver(nelx, nely, volfrac, penal, rmin, ft, gui,
        LBracketBoundaryConditions(nelx, nely, nelx // 3, 2 * nely // 3))
    x_opt  = solver.optimize(x)
    x_opt  = solver.filter_variables(x_opt)
    gui.update(x_opt)

    from PIL import Image
    x_opt.clip(0.0, 1.0)
    x_opt = ((1 - x_opt.reshape(nelx, nely)) * 255).round().astype("uint8").T
    result = Image.fromarray(x_opt)
    result.save("tmp.png")

    # Make sure the plot stays and that the shell remains
    input("Press any key...")

if __name__ == "__main__":
    def run():
        # Default input parameters
        nelx, nely, volfrac, penal, rmin, ft = (sys.argv[1:] +
            [120, 120, 0.2, 6.0, 4, 1][len(sys.argv) - 1:])[:6]
        nelx, nely, ft = map(int, [nelx, nely, ft])
        volfrac, penal, rmin = map(float, [volfrac, penal, rmin])
        main(nelx, nely, volfrac, penal, rmin, ft)
    run()
