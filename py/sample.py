import numpy as np

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import utils as utils


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", ":", "-.", "--"]


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
plt.rcParams['font.size'] = 16

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


# prepare color palette
Ncol, col = utils.set_color_palette_for_color_universal_design()


fig = utils.set_figure(1, 1)
ax = [0]
utils.locate_panels(fig, ax, 1, 1, True, True)


xmin, xmax = 0.0, np.pi
num = 128
xx = np.linspace(xmin, xmax, num)
yy = [0] * num

for ii in range(Ncol):
    yy = np.sin(xx) + 0.2 * ii
    # ax[0].plot(xx, yy, linestyle = ls[ii % Nls], color = col[ii % Ncol])
    ax[0].plot(xx, yy, linestyle = ls[0], color = col[ii % Ncol])

ax[0].set_xlabel(r"$x$")
ax[0].set_ylabel(r"$y$")

ax[0].grid()


fig.savefig("fig/" + "color_sample" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
fig.savefig("fig/" + "color_sample" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.close("all")
