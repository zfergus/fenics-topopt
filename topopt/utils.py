from __future__ import division

import numpy as np


def xy_to_id(x, y, nelx, nely, order="F"):
    if order == "C":
        return (y * (nelx + 1)) + x
    else:
        return (x * (nely + 1)) + y


def id_to_xy(index, nelx, nely, order="F"):
    if order == "C":
        y = index // (nelx + 1)
        x = index % (nelx + 1)
    else:
        x = index // (nely + 1)
        y = index % (nely + 1)
    return x, y


def deleterowcol(A, delrow, delcol):
    """Assumes that matrix is in symmetric csc form !"""
    m = A.shape[0]
    keep = np.delete(np.arange(0, m), delrow)
    A = A[keep, :]
    keep = np.delete(np.arange(0, m), delcol)
    A = A[:, keep]
    return A
