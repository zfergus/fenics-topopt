from __future__ import division

import numpy as np
import scipy


class Filter(object):

    def __init__(self, nelx, nely, rmin):
        """
        Filter: Build (and assemble) the index+data vectors for the coo matrix
        format.
        """
        nfilter = int(nelx * nely * ((2 * (np.ceil(rmin) - 1) + 1)**2))
        iH = np.zeros(nfilter)
        jH = np.zeros(nfilter)
        sH = np.zeros(nfilter)
        cc = 0
        for i in range(nelx):
            for j in range(nely):
                row = i * nely + j
                kk1 = int(np.maximum(i - (np.ceil(rmin) - 1), 0))
                kk2 = int(np.minimum(i + np.ceil(rmin), nelx))
                ll1 = int(np.maximum(j - (np.ceil(rmin) - 1), 0))
                ll2 = int(np.minimum(j + np.ceil(rmin), nely))
                for k in range(kk1, kk2):
                    for l in range(ll1, ll2):
                        col = k * nely + l
                        fac = rmin - np.sqrt(
                            ((i - k) * (i - k) + (j - l) * (j - l)))
                        iH[cc] = row
                        jH[cc] = col
                        sH[cc] = np.maximum(0.0, fac)
                        cc = cc + 1
        # Finalize assembly and convert to csc format
        self.H = scipy.sparse.coo_matrix((sH, (iH, jH)),
            shape=(nelx * nely, nelx * nely)).tocsc()
        self.Hs = self.H.sum(1)

    def filter_variables(self, x, xPhys, ft):
        if ft == 0:
            xPhys[:] = x
        elif ft == 1:
            xPhys[:] = np.asarray(self.H * x[np.newaxis].T / self.Hs)[:, 0]

    def filter_compliance_sensitivities(self, xPhys, dc, ft):
        if ft == 0:
            dc[:] = (np.asarray((self.H * (xPhys * dc))[np.newaxis].T /
                self.Hs)[:, 0] / np.maximum(0.001, xPhys))
        elif ft == 1:
            dc[:] = np.asarray(self.H * (dc[np.newaxis].T / self.Hs))[:, 0]

    def filter_volume_sensitivities(self, _xPhys, dv, ft):
        if ft == 0:
            pass
        elif ft == 1:
            dv[:] = np.asarray(self.H * (dv[np.newaxis].T / self.Hs))[:, 0]
