from __future__ import print_function, division
from fenics import *

from linear_elasticity import sigma


def von_Mises_stress(mesh, displacements):
    """ Computes the stress on a mesh given the displacements. """
    V = FunctionSpace(mesh, "P", 1)
    d = displacements.geometric_dimension()
    # Deviatoric stress
    s = (sigma(displacements) - (1. / 3) * tr(sigma(displacements)) *
        Identity(d))
    von_Mises = sqrt(3. / 2 * inner(s, s))
    von_Mises = project(von_Mises, V)
    return von_Mises
