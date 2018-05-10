import matplotlib.pyplot as plt
from scipy.ndimage import imread
import numpy as np
import sys

from lxml import etree


def mesh_from_img(img):
    nv = (img.shape[0] + 1) * (img.shape[1] + 1)
    nf = img.size * 2

    v_count = 0
    f_count = 0
    V_dict = {}

    V = np.zeros([nv, 2])
    F = np.zeros([nf, 3], dtype=np.int)

    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            val = img[i, j]
            if val == 255.0:
                continue

            v_idx = []
            for v_i in [(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)]:
                if v_i in V_dict:
                    v_idx.append(V_dict[v_i])
                else:
                    V_dict[v_i] = v_count
                    V[v_count, :] = np.array((v_i[1], -v_i[0]))
                    v_count += 1
                    v_idx.append(v_count - 1)

            v1, v2, v3, v4 = v_idx

            F[f_count, :] = np.array([v1, v2, v4])
            F[f_count + 1, :] = np.array([v1, v4, v3])
            f_count += 2

    V = np.resize(V, [v_count, 2])
    F = np.resize(F, [f_count, 3])

    return V, F


def plot_mesh(V, F, plot_f=False):
    plt.scatter(V[:, 0], V[:, 1])
    if not plot_f:
        return
    for i in range(F.shape[0]):
        vf = V[F[i, :], :]
        plt.plot(vf[:, 0], vf[:, 1])


def dump_dolfin_xml(filename, V, F):
    root = etree.Element("dolfin")
    mesh = etree.SubElement(root, "mesh", celltype="triangle", dim="2")
    vertices = etree.SubElement(mesh, "vertices", size="%d" % V.shape[0])
    for i in range(V.shape[0]):
        etree.SubElement(vertices, "vertex", index="%d" % i, x="%f" % V[i, 0],
            y="%f" % V[i, 1], z="0")
    faces = etree.SubElement(mesh, "cells", size="%d" % F.shape[0])
    for i in range(F.shape[0]):
        etree.SubElement(faces, "triangle", index="%d" % i, v0="%d" % F[i, 0],
            v1="%d" % F[i, 1], v2="%d" % F[i, 2])

    with open(filename, 'wb') as f:
        f.write(etree.tostring(root, pretty_print=True))


if __name__ == "__main__":
    img = imread(sys.argv[1])

    plt.figure()
    plt.imshow(img)
    plt.colorbar()

    plt.figure()
    V, F = mesh_from_img(img)
    dump_dolfin_xml("derp.xml", V, F)
    plot_mesh(V, F)
    plt.show()
