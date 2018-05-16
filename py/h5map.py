import h5py
import numpy as np
import math
import time # for performance measurement

import multiprocessing as mp

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm # for logarithmic plot in imshow

import utils as utils


nxpanel, nypanel = 4, 4
filename = "cb17"
tag = ["6", "5", "4", "3", "2", "1", "0.6", "0.3"]
zd  = [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.6, 0.3]


def make_map(idx):
    ii = int(np.floor(idx / nypanel))
    jj =              idx % nypanel

    # pick up an appropriate snapshot
    if (ii & 1) == 0:
        snapshot = "000"
    else:
        snapshot = "040"
    kk = jj * (nxpanel >> 1) + ((nxpanel - 1 - ii) >> 1)
    model = tag[kk] + "kpc"
    input_file = model + "/dat/" + filename + ".split" + snapshot + ".h5"

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

    dSinv = 1.0 / (dx * dz)
    ff *= dSinv

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

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# generate additional color maps
my_cmap = utils.generate_cmap(["darkblue", "deepskyblue", "lime", "yellow", "red", "magenta", "white"])
# cmap0 = utils.generate_cmap(["blue", "cyan", "green", "yellow", "red", "magenta"])
# cmap1 = utils.generate_cmap(["black", "blue", "cyan", "green", "yellow", "red", "magenta", "white"])
# cmap2 = utils.generate_cmap(["blue", "cyan", "green", "yellow", "red", "magenta", "white"])
# cmap3 = utils.generate_cmap(["blue", "green", "yellow", "red", "magenta"])
# cmap4 = utils.generate_cmap(["black", "blue", "green", "yellow", "red", "magenta", "white"])
# cmap5 = utils.generate_cmap(["blue", "green", "yellow", "red", "magenta", "white"])


fig = utils.set_figure(nxpanel, nypanel)
ax = [0] * nxpanel * nypanel
utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

# set plot range
# xmin, xmax = -25.0, 25.0
# zmin, zmax = -25.0, 25.0
# xmin, xmax = -15.0, 15.0
# zmin, zmax = -15.0, 15.0
xmin, xmax = -10.0, 10.0
zmin, zmax = -10.0, 10.0

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
# fmin, fmax = 1.0e+5, 3.0e+7
fmin, fmax = 1.0e+7, 3.1e+9

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
    model = tag[kk] + "kpc"
    input_file = model + "/dat/" + filename + ".split" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    length_unit = h5file["/"].attrs["length_astro_unit_name"]
    time_unit = h5file["/"].attrs["time_astro_unit_name"]
    mass_unit = h5file["/"].attrs["mass_astro_unit_name"]
    current_time = h5file["/"].attrs["time"]

    # close the HDF5 file
    h5file.close()

    # plot the data
    # img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = "jet", aspect = "auto")
    img = ax[idx].imshow(ZZ[idx].T, extent = [xmin, xmax, zmin, zmax], origin='lower', interpolation='none', norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

    # overlay contour
    # ax[idx].contour(XX, YY, ZZ[idx].T, 4, colors = "black", norm = LogNorm(vmin = fmin, vmax = fmax))
    # cont = ax[idx].contour(XX, YY, ZZ[idx].T, levels = [3.1e+5, 1.0e+6, 3.1e+6, 1.0e+7], colors = "black", norm = LogNorm(vmin = fmin, vmax = fmax))
    # cont = ax[idx].contour(XX, YY, ZZ[idx].T, levels = [3.1e+5, 1.0e+6, 3.1e+6, 1.0e+7, 3.1e+7], colors = "black", norm = LogNorm(vmin = fmin, vmax = fmax), linewidths = 0.7)
    cont = ax[idx].contour(XX, YY, ZZ[idx].T, levels = [3.1e+7, 1.0e+8, 3.1e+8, 1.0e+9], colors = "black", norm = LogNorm(vmin = fmin, vmax = fmax), linewidths = 0.7)
    # ax[idx].clabel(cont, inline = 1, fontsize = 10, fmt = ticker.FuncFormatter(utils.scientific))
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
    # ax[idx].set_xticks([-10, -5, 0, 5, 10])
    # ax[idx].set_yticks([-10, -5, 0, 5, 10])
    ax[idx].set_xticks([-9, -6, -3, 0, 3, 6, 9])
    ax[idx].set_yticks([-9, -6, -3, 0, 3, 6, 9])
    ax[idx].tick_params(axis = "both", direction = "in", color = "white", bottom = "on", top = "on", left = "on", right = "on")

    # set label
    if jj == 0:
        ax[idx].set_xlabel(r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
    if ii == 0:
        ax[idx].set_ylabel(r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))

    # set caption
    caption  = "(" + "{:^c}".format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ")"
    caption += " " + r"${:<} = {:.1f}$ {:<}".format("z_\mathrm{d}", zd[kk], length_unit[0].decode('UTF-8'))
    caption += ", $t = {:.0f}$ {:<}".format(current_time[0] / 1000, "Gyr")
    ax[idx].text(xmin + 0.5, zmax - 2.5, caption, fontsize=11, color = "white")


utils.set_shared_xlabel(ax, r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
utils.set_shared_ylabel(ax, r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))

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
