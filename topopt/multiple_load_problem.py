from __future__ import division

import numpy as np

from problem import Problem


class MultipleLoadProblem(Problem):
    def compute_compliance(self, xPhys, dc):
        # Compute compliance and its gradient
        obj = 0.0
        dc[:] = 0.0
        for i in range(self.f.shape[1]):
            ui = self.u[:, i][self.edofMat].reshape(-1, 8)
            self.ce[:] = (ui.dot(self.KE) * ui).sum(1)
            obj += ((self.Emin + xPhys**self.penal *
                (self.Emax - self.Emin)) * self.ce).sum()
            dc[:] += (-self.penal * xPhys**(self.penal - 1.0) *
                (self.Emax - self.Emin)) * self.ce
        dc /= float(self.f.shape[1])
        return obj / float(self.f.shape[1])
