from __future__ import division

import numpy as np
import nlopt

from filter import Filter
from problem import Problem as TopoptProblem
from boundary_conditions import BoundaryConditions


class Solver(object):

    def __init__(self, nelx, nely, volfrac, penal, rmin, ft, gui, bc):
        self.n = nelx * nely
        self.opt = nlopt.opt(nlopt.LD_MMA, self.n)
        self.passive = bc.get_passive_elements()
        self.xPhys = np.ones(self.n)
        if self.passive is not None:
            self.xPhys[self.passive] = 0

        # set bounds
        ub = np.ones(self.n, dtype=float)
        self.opt.set_upper_bounds(ub)
        lb = np.zeros(self.n, dtype=float)
        self.opt.set_lower_bounds(lb)

        # set stopping criteria
        self.opt.set_maxeval(2000)
        self.opt.set_ftol_rel(0.001)

        # set objective and constraint functions
        self.opt.set_min_objective(self.compliance_function)
        self.opt.add_inequality_constraint(self.volume_function, 0)

        # setup filter
        self.ft = ft
        self.filtering = Filter(nelx, nely, rmin)

        # setup problem def
        self.init_problem(nelx, nely, penal, bc)
        self.volfrac = volfrac

        # set GUI callback
        self.init_gui(gui)

    def init_problem(self, nelx, nely, penal, bc):
        self.problem = TopoptProblem(nelx, nely, penal, bc)

    def init_gui(self, gui):
        self.gui = gui
        self.gui.plot_force_arrows(self.problem.f)

    def optimize(self, x):
        self.xPhys = x.copy()
        x = self.opt.optimize(x)
        return x

    def filter_variables(self, x):
        self.filtering.filter_variables(x, self.xPhys, self.ft)
        if self.passive is not None:
            self.xPhys[self.passive] = 0
        return self.xPhys

    def compliance_function_fdiff(self, x, dc):
        obj = self.compliance_function(x, dc)

        x0 = x.copy()
        dc0 = dc.copy()
        dcf = np.zeros(dc.shape)
        for i, v in enumerate(x):
            x = x0.copy()
            x[i] += 1e-6
            o1 = self.compliance_function(x, dc)
            x[i] = x0[i] - 1e-6
            o2 = self.compliance_function(x, dc)
            dcf[i] = (o1 - o2) / (2e-6)
        print("finite differences: {:g}".format(np.linalg.norm(dcf - dc0)))
        dc[:] = dc0

        return obj

    def compliance_function(self, x, dc):
        # Filter design variables
        self.filter_variables(x)

        # Display physical variables
        self.gui.update(self.xPhys)

        # Setup and solve FE problem
        self.problem.compute_displacements(self.xPhys)

        # Objective and sensitivity
        obj = self.problem.compute_compliance(self.xPhys, dc)

        # Sensitivity filtering
        self.filtering.filter_compliance_sensitivities(self.xPhys, dc, self.ft)

        return obj

    def volume_function(self, x, dv):
        # Filter design variables
        self.filter_variables(x)

        # Volume sensitivities
        dv[:] = 1.0

        # Sensitivity filtering
        self.filtering.filter_volume_sensitivities(self.xPhys, dv, self.ft)

        return sum(self.xPhys) - self.volfrac * len(x)
