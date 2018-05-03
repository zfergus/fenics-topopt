# -*- coding: utf-8 -*-
from __future__ import print_function, division

import pdb

import numpy
import scipy.sparse

from utils import xy_to_id, id_to_xy


class VonMisesStressCalculator:
    def __init__(self, nelx, nely, Emin, Emax, penal):
        self.nelx, self.nely, self.Emin, self.Emax, self.penal = (
            nelx, nely, Emin, Emax, penal)
        self.edofMat = self.build_indices()

    @staticmethod
    def B(side):
        """ Precomputed strain-displacement matrix. """
        n = -0.5 / side
        p = 0.5 / side
        return numpy.array([[p, 0, n, 0, n, 0, p, 0],
                            [0, p, 0, p, 0, n, 0, n],
                            [p, p, p, n, n, n, n, p]])

    @staticmethod
    def E(nu):
        """ Precomputed constitutive matrix. """
        return numpy.array([[1, nu, 0],
                            [nu, 1, 0],
                            [0, 0, (1 - nu) / 2.]]) / (1 - nu**2)

    def build_indices(self):
        edofMat = numpy.zeros((8, self.nelx * self.nely), dtype=int)
        for elx in range(self.nelx):
            for ely in range(self.nely):
                el = ely + elx * self.nely # Element index
                n1 = (self.nely + 1) * elx + ely # Left nodes
                n2 = (self.nely + 1) * (elx + 1) + ely # Right nodes
                edofMat[:, el] = numpy.array([2 * n1 + 2, 2 * n1 + 3,
                    2 * n2 + 2, 2 * n2 + 3, 2 * n2, 2 * n2 + 1, 2 * n1,
                    2 * n1 + 1])
        return edofMat

    def penalized_densities(self, x):
        """ Compute the penalized densties. """
        return self.Emin + (self.Emax - self.Emin) * x**self.penal

    def diff_penalized_densities(self, x):
        """ Compute the penalized densties. """
        return (self.Emax - self.Emin) * self.penal * x**(self.penal - 1)

    def calculate_principle_stresses(self, x, u, nu, side=1):
        """
        Calculate the principle stresses in the x, y, and shear directions.
        """
        rho = self.penalized_densities(x)
        EB = self.E(nu).dot(self.B(side))
        stress = sum([EB.dot(u[:, i][self.edofMat]) for i in range(u.shape[1])])
        stress *= rho / float(u.shape[1])
        return numpy.hsplit(stress.T, 3)

    def calculate_stress(self, x, u, nu, side=1):
        """
        Calculate the Von Mises stress given the densities x, displacements u,
        and young modulus nu.
        """
        s11, s22, s12 =  self.calculate_principle_stresses(x, u, nu, side)
        vm_stress = numpy.sqrt(s11**2 - s11 * s22 + s22**2 + 3 * s12**2)
        return vm_stress

    def calculate_diff_stress(self, x, u, nu, side=1):
        """
        Calculate the derivative of the Von Mises stress given the densities x,
        displacements u, and young modulus nu. Optionally, provide the side
        length (default: 1).
        """
        rho = self.penalized_densities(x)
        EB = self.E(nu).dot(self.B(side))
        EBu = sum([EB.dot(u[:, i][self.edofMat]) for i in range(u.shape[1])])
        s11, s22, s12 = numpy.hsplit((EBu * rho / float(u.shape[1])).T, 3)
        drho = self.diff_penalized_densities(x)
        ds11, ds22, ds12 = numpy.hsplit(
            ((1 - rho) * drho * EBu / float(u.shape[1])).T, 3)
        vm_stress = numpy.sqrt(s11**2 - s11 * s22 + s22**2 + 3 * s12**2)
        if abs(vm_stress).sum() > 1e-8:
            dvm_stress = (0.5 * (1. / vm_stress) * (2 * s11 * ds11 -
                ds11 * s22 - s11 * ds22 + 2 * s22 * ds22 + 6 * s12 * ds12))
            return dvm_stress
        return 0

    def calculate_fdiff_stress(self, x, u, nu, side=1, dx=1e-6):
        """
        Calculate the derivative of the Von Mises stress using finite
        differences given the densities x, displacements u, and young modulus
        nu. Optionally, provide the side length (default: 1) and delta x
        (default: 1e-6).
        """
        ds = self.calculate_diff_stress(x, u, nu, side)
        dsf = numpy.zeros(x.shape)
        x = numpy.expand_dims(x, -1)
        for i in range(x.shape[0]):
            delta = scipy.sparse.coo_matrix(([dx], [[i], [0]]), shape=x.shape)
            s1 = self.calculate_stress((x + delta.A).squeeze(), u, nu, side)
            s2 = self.calculate_stress((x - delta.A).squeeze(), u, nu, side)
            dsf[i] = ((s1 - s2) / (2. * dx))[i]
        print("finite differences: {:g}".format(numpy.linalg.norm(dsf - ds)))
        return dsf
