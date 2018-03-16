import matplotlib.pyplot as plt


# set figure size and its aspect ratio
def set_figure(nxpanel, nypanel):
    Lx = 10
    if nxpanel > 2 * nypanel:
        Lx *= 2
    Ly = (Lx / nxpanel) * nypanel
    fig = plt.figure(figsize = (Lx, Ly))

    return fig


def locate_panels(fig, ax, nx, ny, share_xaxis, share_yaxis):
    margin = 0.12
    if (share_xaxis == False) or (share_yaxis == False):
        margin = 0.15

    xmin, xmax = margin, 1.0 - margin
    ymin, ymax = margin, 1.0 - margin
    xbin = (xmax - xmin) / nx
    ybin = (ymax - ymin) / ny
    xmargin, ymargin = 0, 0

    if share_yaxis == False:
        xmin = 0.0
        xbin = 1.0 / nx
        xmargin = xbin * margin

    if share_xaxis == False:
        ymin = 0.0
        ybin = 1.0 / ny
        ymargin = ybin * margin

    for ii in range(nx):
        xl = xmin + ii * xbin + xmargin

        for jj in range(ny):
            yl = ymin + jj * ybin + ymargin
            kk = ii * ny + jj
            ax[kk] = fig.add_axes((xl, yl, xbin - 2 * xmargin, ybin - 2 * ymargin))
            ax[kk].tick_params(axis = "both", direction = "in", bottom = True, top = True, left = True, right = True)

            if share_xaxis == True:
                ax[kk].tick_params(labelbottom = False)
                if jj == 0:
                    ax[kk].tick_params(labelbottom = True)

            if share_yaxis == True:
                ax[kk].tick_params(labelleft = False)
                if ii == 0:
                    ax[kk].tick_params(labelleft = True)


def set_shared_xlabel(ax, xlabel):
    fig = ax[-1].figure
    fig.canvas.draw()

    # get the corner for all plots
    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax:
        at.set_xlabel("") # remove existing xlabels

        xl, xr = at.get_position().x0, at.get_position().x1

        bboxes, _ = at.xaxis.get_ticklabel_extents(fig.canvas.renderer)
        bboxes = bboxes.inverse_transformed(fig.transFigure)
        yb, yt = bboxes.y0, bboxes.y1

        if x0 > xl:
            x0 = xl
        if x1 < xr:
            x1 = xr
        if y0 > yb:
            y0 = yb
        if y1 < yt:
            y1 = yt

    # set position of label
    ax[-1].set_xlabel(xlabel)
    # ax[-1].xaxis.set_label_coords((x0 + x1) / 2, (y0 + y1) / 2, transform = fig.transFigure)
    ax[-1].xaxis.set_label_coords((x0 + x1) / 2, y0, transform = fig.transFigure)


def set_shared_ylabel(ax, ylabel):
    fig = ax[-1].figure
    fig.canvas.draw()

    # get the corner for all plots
    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax:
        at.set_ylabel("") # remove existing ylabels

        yb, yt = at.get_position().y0, at.get_position().y1

        bboxes, _ = at.yaxis.get_ticklabel_extents(fig.canvas.renderer)
        bboxes = bboxes.inverse_transformed(fig.transFigure)
        xl, xr = bboxes.x0, bboxes.x1

        if x0 > xl:
            x0 = xl
        if x1 < xr:
            x1 = xr
        if y0 > yb:
            y0 = yb
        if y1 < yt:
            y1 = yt

    # set position of label
    ax[-1].set_ylabel(ylabel)
    ax[-1].yaxis.set_label_coords((x0 + x1) / 2, (y0 + y1) / 2, transform = fig.transFigure)


# for scientific notation in coutour plot
# based on http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib
import matplotlib.ticker as ticker
def scientific(x, pos):
    a, b = "{:.1e}".format(x).split("e")
    b = int(b)
    return r"${} \times 10^{{{}}}$".format(a, b)


# generate color map
import sys
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
def generate_cmap(cols):
    vals = range(len(cols))
    vmax = np.ceil(np.max(vals))
    color_list = []
    for vv, cc in zip(vals, cols):
        color_list.append((vv / vmax, cc))
    return LinearSegmentedColormap.from_list("custom_cmap", color_list)

def generate_cmap_comp(cols, loci):
    if len(cols) != len(loci):
        print("colors and position list need the same number of elements.")
        sys.exit(0)
    vals = range(len(cols))
    vmax = np.ceil(np.max(vals))
    color_list = []
    for pp, cc in zip(loci, cols):
        color_list.append((pp, cc))
    return LinearSegmentedColormap.from_list("custom_cmap", color_list)
