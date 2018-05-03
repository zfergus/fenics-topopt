from __future__ import division

import numpy as np

import pdb

import scipy.sparse
from scipy.sparse import coo_matrix
import cvxopt
import cvxopt.cholmod

from utils import deleterowcol


class Problem(object):

    @staticmethod
    def lk(E=1.):
        """element stiffness matrix"""
        nu = 0.3
        k = np.array([0.5 - nu / 6., 0.125 + nu / 8., -0.25 - nu / 12.,
            -0.125 + 0.375 * nu, -0.25 + nu / 12., -0.125 - nu / 8., nu / 6.,
            0.125 - 0.375 * nu])
        KE = E / (1 - nu**2) * np.array([
            [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
            [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
            [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
            [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
            [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
            [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
            [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
            [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]])
        return KE

    def build_indices(self, nelx, nely):
        """ FE: Build the index vectors for the for coo matrix format. """
        self.KE = self.lk()
        self.edofMat = np.zeros((nelx * nely, 8), dtype=int)
        for elx in range(nelx):
            for ely in range(nely):
                el = ely + elx * nely
                n1 = (nely + 1) * elx + ely
                n2 = (nely + 1) * (elx + 1) + ely
                self.edofMat[el, :] = np.array([2 * n1 + 2, 2 * n1 + 3,
                    2 * n2 + 2, 2 * n2 + 3, 2 * n2, 2 * n2 + 1, 2 * n1,
                    2 * n1 + 1])
        # Construct the index pointers for the coo format
        self.iK = np.kron(self.edofMat, np.ones((8, 1))).flatten()
        self.jK = np.kron(self.edofMat, np.ones((1, 8))).flatten()

    def __init__(self, nelx, nely, penal, bc):
        # Problem size
        self.nelx = nelx
        self.nely = nely

        # Max and min stiffness
        self.Emin = 1e-9
        self.Emax = 1.0

        # SIMP penalty
        self.penal = penal

        # dofs:
        self.ndof = 2 * (nelx + 1) * (nely + 1)

        # FE: Build the index vectors for the for coo matrix format.
        self.build_indices(nelx, nely)

        # BC's and support (half MBB-beam)
        dofs = np.arange(2 * (nelx + 1) * (nely + 1))
        self.fixed = bc.get_fixed_nodes()
        self.free = np.setdiff1d(dofs, self.fixed)

        # Solution and RHS vectors
        self.f = bc.get_forces()
        self.u = np.zeros(self.f.shape)

        # Per element compliance
        self.ce = np.zeros(nely * nelx)

    def compute_displacements(self, xPhys):
        # Setup and solve FE problem
        sK = ((self.KE.flatten()[np.newaxis]).T * (
            self.Emin + (xPhys)**self.penal *
            (self.Emax - self.Emin))).flatten(order='F')
        K = scipy.sparse.coo_matrix((sK, (self.iK, self.jK)),
            shape=(self.ndof, self.ndof)).tocsc()
        # Remove constrained dofs from matrix and convert to coo
        K = deleterowcol(K, self.fixed, self.fixed).tocoo()
        # Solve system
        K1 = cvxopt.spmatrix(K.data, K.row.astype(np.int), K.col.astype(np.int))
        B = cvxopt.matrix(self.f[self.free, :])
        cvxopt.cholmod.linsolve(K1, B)
        self.u[self.free, :] = np.array(B)[:, :]

    def compute_compliance(self, xPhys, dc):
        # Compute compliance and its gradient
        u = self.u[:, 0][self.edofMat].reshape(-1, 8)
        self.ce[:] = (u.dot(self.KE) * u).sum(1)
        obj = ((self.Emin + xPhys**self.penal *
            (self.Emax - self.Emin)) * self.ce).sum()
        dc[:] = (-self.penal * xPhys**(self.penal - 1.0) *
            (self.Emax - self.Emin)) * self.ce
        return obj
