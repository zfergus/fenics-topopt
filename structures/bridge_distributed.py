#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import pdb

import numpy as np

from topopt.gui import GUI
import topopt.boundary_conditions
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

# BoundaryConditions = topopt.boundary_conditions.BoundaryConditions


class BoundaryConditions(topopt.boundary_conditions.BoundaryConditions):
    def get_forces(self):
        # Return the force vector for the problem
        topx_to_id = np.vectorize(
            lambda x: xy_to_id(x, 0, self.nelx, self.nely))
        topx = 2 * topx_to_id(np.arange((self.nelx + 1) // 2)) + 1
        f = np.zeros((2 * (self.nelx + 1) * (self.nely + 1), 1))
        f[topx, 0] = -100
        return f


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
        BoundaryConditions(nelx, nely))
    x_opt  = solver.optimize(x)
    x_opt  = solver.filter_variables(x_opt)
    gui.update(x_opt)

    from PIL import Image
    x_opt.clip(0.0, 1.0)
    x_opt = ((1 - x_opt.reshape(nelx, nely)) * 255).round().astype("uint8").T
    x_opt_sym = x_opt[:, ::-1]
    result = Image.fromarray(np.hstack([x_opt_sym, x_opt]))
    result.save("tmp.png")

    # Make sure the plot stays and that the shell remains
    input("Press any key...")

if __name__ == "__main__":
    def run():
        # Default input parameters
        nelx, nely, volfrac, penal, rmin, ft = (sys.argv[1:] +
            [300, 100, 0.4, 3.0, 5.4, 1][len(sys.argv) - 1:])[:6]
        nelx, nely, ft = map(int, [nelx, nely, ft])
        volfrac, penal, rmin = map(float, [volfrac, penal, rmin])
        main(nelx, nely, volfrac, penal, rmin, ft)
    run()
