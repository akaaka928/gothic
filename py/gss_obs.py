import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm # for logarithmic plot in imshow

import os.path as path
import multiprocessing as mp

import utils as utils
import m31 as m31


nxpanel, nypanel = 1, 1

filename = "m16disk"
init = 0
last = 399


xmin, xmax = -4.0, 8.0
ymin, ymax = -8.0, 4.0
fmin, fmax = 1.0e+4, 1.0e+7

nx, ny = 256, 256
dx, dy = (xmax - xmin) / nx, (ymax - ymin) / ny

# settings for point spread function with gaussian
smooth = 1
# smooth = 1.5
# smooth = 2
# contain = 1
# contain = 2
contain = 3
sigma = smooth * (dx + dy) / 2
PSFinv = 1 / (math.sqrt(2) * sigma)
Nsmooth = int(np.ceil(contain * smooth))


def draw_figure(fileid, rot):
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".snp" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    time_unit = h5file["/axis labels/"].attrs["time_astro_unit_name"]
    time = h5file["/snp/"].attrs["time"]

    xx = np.linspace(xmin, xmax, nx)
    yy = np.linspace(ymin, ymax, ny)

    dSxyinv = 1.0 / (dx * dy)

    # read particle position and mass
    folder = "snp/"
    position = h5file[folder + "position"].value
    mass = h5file[folder + "mass"].value

    # data preparation
    px = position[:, 0]
    py = position[:, 1]
    pz = position[:, 2]

    xy_map = np.zeros((nx, ny))

    for jj in range(mass.size):
        # coordinate transformation (disk orthogonal frame ==>> observed frame)
        pos = np.array([px[jj], py[jj], pz[jj]])
        obs = np.dot(rot, pos)

        # coordinate transformation (M31 standard coordinate)
        xi, eta, D = m31.standard_coordinate(obs[0], obs[1], obs[2])

        xj = int(np.rint((xi  - xmin) / dx - 0.5))
        yj = int(np.rint((eta - ymin) / dy - 0.5))

        # apply smoothing by point spreading function
        xp = np.linspace(xmin + (xj - Nsmooth) * dx, xmin + (xj + Nsmooth + 1) * dx, 2 * Nsmooth + 2)
        yp = np.linspace(ymin + (yj - Nsmooth) * dy, ymin + (yj + Nsmooth + 1) * dy, 2 * Nsmooth + 2)

        erfx = [0] * (2 * Nsmooth + 2)
        erfy = [0] * (2 * Nsmooth + 2)
        for kk in range(2 * Nsmooth + 2):
            erfx[kk] = math.erf((xp[kk] -  xi) * PSFinv)
            erfy[kk] = math.erf((yp[kk] - eta) * PSFinv)

        psf_x = [0] * (2 * Nsmooth + 1)
        psf_y = [0] * (2 * Nsmooth + 1)
        for kk in range(2 * Nsmooth + 1):
            psf_x[kk] = 0.5 * (erfx[kk + 1] - erfx[kk])
            psf_y[kk] = 0.5 * (erfy[kk + 1] - erfy[kk])

        for kk in range(max(0, xj - Nsmooth), min(nx - 1, xj + Nsmooth)):
            for ll in range(max(0, yj - Nsmooth), min(ny - 1, yj + Nsmooth)):
                xy_map[kk][ll] += mass[jj] * psf_x[kk - xj - Nsmooth] * psf_y[ll - yj - Nsmooth]


    # close the HDF5 file
    h5file.close()

    xy_map *= dSxyinv
    xy_map = np.maximum(xy_map, fmin)
    xy_map = np.minimum(xy_map, fmax)

    # plot the data
    img = ax[0].imshow(xy_map.T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

    ax[0].set_xlim([xmin, xmax])
    ax[0].set_ylim([ymin, ymax])

    ax[0].spines[  "top"].set_color("white")
    ax[0].spines[ "left"].set_color("white")
    ax[0].spines["right"].set_color("white")

    ax[0].set_xlabel(r"$\xi$ ({:<})".format(length_unit[0].decode('UTF-8')))

    ax[0].set_ylabel(r"$y$ ({:<})".format(length_unit[0].decode('UTF-8')))


    ax[0].set_xticks([-2, 0, 2, 4, 6])
    ax[0].set_yticks([-8, -6, -4, -2, 0, 2])
    ax[0].tick_params(axis = "both", direction = "in", color = "white", bottom = "on", top = "on", left = "on", right = "on")


    # add colorbar
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


    # add current time
    # fig.suptitle(r"$t = {:.2f}$ {:<}".format(time[0] / 1000, "Gyr"))
    fig.suptitle(r"$t = {:.0f}$ {:<}".format(time[0], "Myr"))
    # fig.suptitle(r"$t = {:.3f}$ {:<}".format(time[0] / 1000, "Gyr"))


    # plt.show()
    figname = "fig/" + filename + "_obs" + snapshot
    plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")


def wrapper(argv):
    return draw_figure(*argv)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 14


# set rotation matrix between observed coordinate and disk orthogonal coordinate
rot, inv = m31.Euler_angle()

# set reference points of shells and GSS
Neast, Eshell_xi, Eshell_eta, Nwest, Wshell_xi, Wshell_eta = m31.Fardal2007_shell()
Nfield, field_xi, field_eta = m31.GSS_obs_field()


my_cmap = utils.generate_cmap(["darkblue", "deepskyblue", "lime", "yellow", "red", "magenta", "white"])

cores = mp.cpu_count()
pool = mp.Pool(cores)
args = [(ii, inv) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
