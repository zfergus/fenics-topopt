from __future__ import division

from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

from utils import xy_to_id, id_to_xy


class GUI(object):
    """
    Graphics user interface for drawing the outputs a topology optimization
    problem
    """
    def __init__(self, nelx, nely, title=""):
        """Initialize plot and plot the initial design"""
        plt.ion()  # Ensure that redrawing is possible
        self.fig, self.ax = plt.subplots()
        self.im = self.ax.imshow(-np.zeros((nelx, nely)).T, cmap='gray',
            interpolation='none', norm=colors.Normalize(vmin=-1, vmax=0))
        plt.xlabel(title)
        # self.fig.tight_layout()
        self.fig.show()
        self.nelx, self.nely = nelx, nely

    def plot_force_arrows(self, f):
        """Add arrows to the plot for each force."""
        arrowprops = {"arrowstyle": "->", "connectionstyle": "arc3",
            "lw": "2", "color": 0}
        cmap = plt.cm.get_cmap("hsv", f.shape[1] + 1)
        for load_i in range(f.shape[1]):
            nz = np.nonzero(f[:, load_i])
            arrowprops["color"] = cmap(load_i)
            for i in range(nz[0].shape[0]):
                x, y = id_to_xy(nz[0][i] // 2, self.nelx, self.nely)
                x = max(min(x, self.nelx - 1), 0)
                y = max(min(y, self.nely - 1), 0)
                z = int(nz[0][i] % 2)
                mag = -50 * f[nz[0][i], load_i]
                self.ax.annotate("", xy=(x, y), xycoords="data",
                    xytext = (0 if z else mag, mag if z else 0),
                    textcoords="offset points", arrowprops=arrowprops)

    def update(self, xPhys, title=None):
        """Plot to screen"""
        self.im.set_array(-xPhys.reshape((self.nelx, self.nely)).T)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        if title is not None:
            plt.title(title)
        plt.pause(0.01)
