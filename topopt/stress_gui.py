from __future__ import division

from matplotlib import colors
import matplotlib.cm as colormaps
import matplotlib.pyplot as plt
import numpy as np

from gui import GUI
from utils import xy_to_id, id_to_xy


class StressGUI(GUI):
    """
    Graphics user interface for drawing the outputs a topology optimization
    problem
    """
    def __init__(self, nelx, nely, stress_calculator, nu, title=""):
        """Initialize plot and plot the initial design"""
        super(StressGUI, self).__init__(nelx, nely, title)
        self.stress_im = self.ax.imshow(
            np.swapaxes(np.zeros((nelx, nely, 4)), 0, 1),
            norm=colors.Normalize(vmin=0, vmax=1), cmap='jet')
        self.fig.colorbar(self.stress_im)
        self.stress_calculator = stress_calculator
        self.nu = nu
        self.myColorMap = colormaps.ScalarMappable(
            norm=colors.Normalize(vmin=0, vmax=1), cmap=colormaps.jet)

    def update(self, xPhys, u, title=None):
        """Plot to screen"""
        self.im.set_array(-xPhys.reshape((self.nelx, self.nely)).T)
        stress = self.stress_calculator.calculate_stress(xPhys, u, self.nu)
        # self.stress_calculator.calculate_fdiff_stress(xPhys, u, self.nu)
        self.myColorMap.set_norm(colors.Normalize(vmin=0, vmax=max(stress)))
        stress_rgba = self.myColorMap.to_rgba(stress)
        stress_rgba[:, :, 3] = xPhys.reshape(-1, 1)
        self.stress_im.set_array(np.swapaxes(
            stress_rgba.reshape((self.nelx, self.nely, 4)), 0, 1))
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        if title is not None:
            plt.title(title)
        else:
            plt.xlabel("Max stress = {:.2f}".format(max(stress)[0]))
        plt.pause(0.01)
