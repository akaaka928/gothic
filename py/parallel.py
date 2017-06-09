import h5py
import numpy as np
import math
import time # for performance measurement

import multiprocessing as mp
# import threading

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm # for logarithmic plot in imshow


# generate color map
import sys
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


# for scientific notation in coutour plot
# based on http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib
import matplotlib.ticker as ticker
def scientific(x, pos):
    a, b = "{:.1e}".format(x).split("e")
    b = int(b)
    return r"${} \times 10^{{{}}}$".format(a, b)


def locate_panels(ax, nx, ny, share_xaxis, share_yaxis):
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

            if share_xaxis == True:
                ax[kk].tick_params(labelbottom = "off")
                if jj == 0:
                    ax[kk].tick_params(labelbottom = "on")

            if share_yaxis == True:
                ax[kk].tick_params(labelleft = "off")
                if ii == 0:
                    ax[kk].tick_params(labelleft = "on")


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
    ax[-1].xaxis.set_label_coords((x0 + x1) / 2, (y0 + y1) / 2, transform = fig.transFigure)


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


def make_map(idx):
    ii = int(np.floor(idx / nypanel))
    jj =              idx % nypanel

    # pick up an appropriate snapshot
    if (ii & 1) == 0:
        snapshot = "000"
    else:
        snapshot = "040"
    kk = jj * (nxpanel >> 1) + ((nxpanel - 1 - ii) >> 1)
    param = ((nxpanel >> 1) * nypanel - 1 - kk) + 1
    model = str(param) + "kpc"
    input_file = model + "/dat/ltg.split" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read particle position and mass
    folder = "data2/"
    position = h5file[folder + "position"].value
    mass = h5file[folder + "mass"].value

    # close the HDF5 file
    h5file.close()

    # data preparation
    px = position[:, 0]
    py = position[:, 1]
    pz = position[:, 2]

    ff = np.zeros((nx, ny))

    for ll in range(mass.size):
        xi = int(np.rint((px[ll] - xmin) / dx - 0.5))
        zj = int(np.rint((pz[ll] - zmin) / dz - 0.5))

        # with smoothing by point spreading function (with performance optimization)
        xp = np.linspace(xmin + (xi - Nsmooth) * dx, xmin + (xi + Nsmooth + 1) * dx, 2 * Nsmooth + 2)
        zp = np.linspace(zmin + (zj - Nsmooth) * dz, zmin + (zj + Nsmooth + 1) * dz, 2 * Nsmooth + 2)
        erfx = [0] * (2 * Nsmooth + 2)
        erfy = [0] * (2 * Nsmooth + 2)
        for mm in range(2 * Nsmooth + 2):
            erfx[mm] = math.erf((xp[mm] - px[ll]) * PSFinv)
            erfy[mm] = math.erf((zp[mm] - pz[ll]) * PSFinv)
        psf_x = [0] * (2 * Nsmooth + 1)
        psf_y = [0] * (2 * Nsmooth + 1)
        for mm in range(2 * Nsmooth + 1):
            psf_x[mm] = 0.5 * (erfx[mm + 1] - erfx[mm])
            psf_y[mm] = 0.5 * (erfy[mm + 1] - erfy[mm])

        for mm in range(max(0, xi - Nsmooth), min(nx - 1, xi + Nsmooth)):
            for nn in range(max(0, zj - Nsmooth), min(ny - 1, zj + Nsmooth)):
                ff[mm][nn] += mass[ll] * psf_x[mm - xi - Nsmooth] * psf_y[nn - zj - Nsmooth]

    ff = np.maximum(ff, fmin)
    ff = np.minimum(ff, fmax)

    return ff


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 14

# generate additional color maps
cmap0 = generate_cmap(["blue", "cyan", "green", "yellow", "red", "magenta"])
cmap1 = generate_cmap(["black", "blue", "cyan", "green", "yellow", "red", "magenta", "white"])
cmap2 = generate_cmap(["blue", "cyan", "green", "yellow", "red", "magenta", "white"])
cmap3 = generate_cmap(["blue", "green", "yellow", "red", "magenta"])
cmap4 = generate_cmap(["black", "blue", "green", "yellow", "red", "magenta", "white"])
cmap5 = generate_cmap(["blue", "green", "yellow", "red", "magenta", "white"])

# set number of panels
nxpanel, nypanel = 4, 4
# nxpanel, nypanel = 2, 2
ax = [0] * nxpanel * nypanel

# set figure size and its aspect ratio
Lx = 10
if nxpanel > 2 * nypanel:
    Lx *= 2
Ly = (Lx / nxpanel) * nypanel
fig = plt.figure(figsize = (Lx, Ly))

# set location of panels
locate_panels(ax, nxpanel, nypanel, True, True)

# set plot range
# xmin, xmax = -25.0, 25.0
# zmin, zmax = -25.0, 25.0
xmin, xmax = -15.0, 15.0
zmin, zmax = -15.0, 15.0

# nx, ny = 32, 32
# nx, ny = 64, 64
# nx, ny = 128, 128
nx, ny = 256, 256
xx = np.linspace(xmin, xmax, nx)
yy = np.linspace(zmin, zmax, ny)
dx, dz = (xmax - xmin) / nx, (zmax - zmin) / ny
XX, YY = np.meshgrid(xx, yy)

# fmin, fmax = 1.1e+5, 1.0e+9
# fmin, fmax = 1.0e+5, 1.0e+8
fmin, fmax = 1.0e+5, 3.0e+7

# settings for point spread function with gaussian
smooth = 1
# smooth = 1.5
# smooth = 2
# contain = 1
# contain = 2
contain = 3
sigma = smooth * (dx + dz) / 2
PSFinv = 1 / (math.sqrt(2) * sigma)
Nsmooth = int(np.ceil(contain * smooth))

# start benchmark
start = time.time()

cores = mp.cpu_count()
cores = int(np.ceil(cores / 2))

# ZZ = np.zeros((nxpanel * nypanel, nx, ny))
# for idx in range(nxpanel * nypanel):
#     make_map(idx)
pool = mp.Pool(cores)
ZZ = pool.map(make_map, range(nxpanel * nypanel))

for idx in range(nxpanel * nypanel):
    ii = int(np.floor(idx / nypanel))
    jj =              idx % nypanel

    # pick up an appropriate snapshot
    if (ii & 1) == 0:
        snapshot = "000"
    else:
        snapshot = "040"
    kk = jj * (nxpanel >> 1) + ((nxpanel - 1 - ii) >> 1)
    param = ((nxpanel >> 1) * nypanel - 1 - kk) + 1
    model = str(param) + "kpc"
    input_file = model + "/dat/ltg.split" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    length_unit = h5file["/"].attrs["length_astro_unit_name"]
    time_unit = h5file["/"].attrs["time_astro_unit_name"]
    mass_unit = h5file["/"].attrs["mass_astro_unit_name"]
    current_time = h5file["/"].attrs["time"]

    # # confirmation of attributes
    # print("{:.1f} {:<}".format(current_time[0], time_unit[0].decode('UTF-8')))
    # print("unit of the length is {:<}".format(length_unit[0].decode('UTF-8')))
    # print("unit of the mass is {:<}".format(mass_unit[0].decode('UTF-8')))

    # close the HDF5 file
    h5file.close()

    # plot the data
    # img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', vmin = fmin, vmax = fmax, cmap = "jet")
    # img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = "jet")
    # img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = cmap0)
    # img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = cmap1)
    # img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = cmap2)
    img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = cmap3)
    # img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = cmap4)
    # img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = cmap5)

    # overlay contour
    # ax[idx].contour(XX, YY, ZZ[idx].T, 4, colors = "black", norm = LogNorm(vmin = fmin, vmax = fmax))
    # cont = ax[idx].contour(XX, YY, ZZ[idx].T, levels = [3.1e+5, 1.0e+6, 3.1e+6, 1.0e+7], colors = "black", norm = LogNorm(vmin = fmin, vmax = fmax))
    # cont = ax[idx].contour(XX, YY, ZZ[idx].T, levels = [3.1e+5, 1.0e+6, 3.1e+6, 1.0e+7, 3.1e+7], colors = "black", norm = LogNorm(vmin = fmin, vmax = fmax), linewidths = 0.7)
    cont = ax[idx].contour(XX, YY, ZZ[idx].T, levels = [3.1e+5, 1.0e+6, 3.1e+6, 1.0e+7], colors = "black", norm = LogNorm(vmin = fmin, vmax = fmax), linewidths = 0.7)
    # ax[idx].clabel(cont, inline = 1, fontsize = 10, fmt = ticker.FuncFormatter(scientific))
    # ax[idx].clabel(cont, inline = 1, fontsize = 10, fmt = "%.1e")

    # # plot color bar
    # if ii == nxpanel - 1:
    #     colorbar_ax = fig.add_axes([ax[idx].get_position().x1, ax[idx].get_position().y0, 0.1 / nxpanel, ax[idx].get_position().y1 - ax[idx].get_position().y0])
    #     cbar = fig.colorbar(img, cax = colorbar_ax, label = r"$\Sigma$ ($M_\odot$ kpc$^{-2}$)")
    #     cbar.solids.set_edgecolor("face")

    # set plot range
    ax[idx].set_xlim([xmin, xmax])
    ax[idx].set_ylim([zmin, zmax])

    # set axis color
    if ii != 0:
        ax[idx].spines["left"].set_color("white")
    if ii != nxpanel - 1:
        ax[idx].spines["right"].set_color("white")
    if jj != 0:
        ax[idx].spines["bottom"].set_color("white")
    if jj != nypanel - 1:
        ax[idx].spines["top"].set_color("white")

    # set ticks
    # ax[idx].set_xticks([-20, -10, 0, 10, 20])
    # ax[idx].set_yticks([-20, -10, 0, 10, 20])
    ax[idx].set_xticks([-10, -5, 0, 5, 10])
    ax[idx].set_yticks([-10, -5, 0, 5, 10])
    ax[idx].tick_params(axis = "both", direction = "in", color = "white", bottom = "on", top = "on", left = "on", right = "on")

    # set label
    if jj == 0:
        ax[idx].set_xlabel(r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
    if ii == 0:
        ax[idx].set_ylabel(r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))

    # set caption
    caption  = "(" + "{:^c}".format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ")"
    caption += " " + r"${:<} = {:.0f}$ {:<}".format("z_\mathrm{d}", param, length_unit[0].decode('UTF-8'))
    caption += ", $t = {:.0f}$ {:<}".format(current_time[0] / 1000, "Gyr")
    ax[idx].text(xmin + 1, zmax - 4, caption, fontsize=12, color = "white")


set_shared_xlabel(ax, r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
set_shared_ylabel(ax, r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))

x0, x1 = 1, 0
y0, y1 = 1, 0
for at in ax:
    xl, xr = at.get_position().x0, at.get_position().x1
    yb, yt = at.get_position().y0, at.get_position().y1

    if x0 > xl:
        x0 = xl
    if x1 < xr:
        x1 = xr
    if y0 > yb:
        y0 = yb
    if y1 < yt:
        y1 = yt
colorbar_ax = fig.add_axes([x1, y0, 0.1 / nxpanel, y1 - y0])
cbar = fig.colorbar(img, cax = colorbar_ax, label = r"$\Sigma$ ($M_\odot$ kpc$^{-2}$)")
cbar.solids.set_edgecolor("face")

# finish benchmark
print("Nsmooth = {:d}, elapsed time is {:e} sec by {:d} CPU cores".format(Nsmooth, time.time() - start, cores))


# plt.show()
plt.savefig("map.png", format = "png", dpi = 300, bbox_inches = "tight")
plt.savefig("map.pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.savefig("map.svg", format = "svg", dpi = 300, bbox_inches = "tight")
