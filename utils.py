from __future__ import print_function, division
from fenics import *
import mshr


def scale_mesh(mesh, width, height):
    """ Scale coordinates of the mesh to be (0, 0) to (width x height). """
    coords = mesh.coordinates()
    x = coords[:, 0]
    y = coords[:, 1]
    x[:] = (x - min(x)) / (max(x) - min(x)) * width
    y[:] = (y - min(y)) / (max(y) - min(y)) * height


def generate_tower_mesh(quality=20):
    otri = mshr.Polygon([Point(0, 0), Point(1, 0), Point(0.5, 1)])
    itri = mshr.Polygon([Point(0.2, 0.0), Point(0.8, 0.0), Point(0.5, 0.75)])
    domain = otri - itri
    return mshr.generate_mesh(domain, quality)


def generate_square_mesh(diagonal="left"):
    return RectangleMesh(Point(0, 0), Point(1, 1), 1, 1, diagonal)
