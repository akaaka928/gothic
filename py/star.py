import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm # for logarithmic plot in imshow

import os.path as path
import multiprocessing as mp

import utils as utils


nxpanel, nypanel = 1, 2

# filename = "cb17"
filename = "m31"
Ndark = 1
init = 0
last = 47
# last = 0
# xmin, xmax = -10.0, 10.0
# ymin, ymax = -10.0, 10.0
# zmin, zmax = -10.0, 10.0
xmin, xmax = -15.0, 15.0
ymin, ymax = -15.0, 15.0
zmin, zmax = -15.0, 15.0
# fmin, fmax = 1.0e+7, 3.1e+9
fmin, fmax = 1.0e+7, 1.0e+10

# nx, ny, nz = 128, 128, 128
nx, ny, nz = 256, 256, 256
# nx, ny, nz = 512, 512, 512

dx, dy, dz = (xmax - xmin) / nx, (ymax - ymin) / ny, (zmax - zmin) / nz

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


def draw_figure(fileid, Nkind):
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".split" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    length_unit = h5file["/"].attrs["length_astro_unit_name"]
    time_unit = h5file["/"].attrs["time_astro_unit_name"]
    mass_unit = h5file["/"].attrs["mass_astro_unit_name"]
    time = h5file["/"].attrs["time"]

    xx = np.linspace(xmin, xmax, nx)
    yy = np.linspace(ymin, ymax, ny)
    zz = np.linspace(zmin, zmax, nz)

    dSxyinv = 1.0 / (dx * dy)
    dSxzinv = 1.0 / (dx * dz)

    xy_sum = np.zeros((nx, ny))
    xz_sum = np.zeros((nx, nz))


    for ii in range(Ndark, Nkind):
        # read particle position and mass
        folder = "data" + str(ii) + "/"
        position = h5file[folder + "position"].value
        mass = h5file[folder + "mass"].value

        # data preparation
        px = position[:, 0]
        py = position[:, 1]
        pz = position[:, 2]

        xy_map = np.zeros((nx, ny))
        xz_map = np.zeros((nx, nz))

        for jj in range(mass.size):
            xj = int(np.rint((px[jj] - xmin) / dx - 0.5))
            yj = int(np.rint((py[jj] - ymin) / dy - 0.5))
            zj = int(np.rint((pz[jj] - zmin) / dz - 0.5))

            # apply smoothing by point spreading function
            xp = np.linspace(xmin + (xj - Nsmooth) * dx, xmin + (xj + Nsmooth + 1) * dx, 2 * Nsmooth + 2)
            yp = np.linspace(ymin + (yj - Nsmooth) * dy, ymin + (yj + Nsmooth + 1) * dy, 2 * Nsmooth + 2)
            zp = np.linspace(zmin + (zj - Nsmooth) * dz, zmin + (zj + Nsmooth + 1) * dz, 2 * Nsmooth + 2)

            erfx = [0] * (2 * Nsmooth + 2)
            erfy = [0] * (2 * Nsmooth + 2)
            erfz = [0] * (2 * Nsmooth + 2)
            for kk in range(2 * Nsmooth + 2):
                erfx[kk] = math.erf((xp[kk] - px[jj]) * PSFinv)
                erfy[kk] = math.erf((yp[kk] - py[jj]) * PSFinv)
                erfz[kk] = math.erf((zp[kk] - pz[jj]) * PSFinv)

            psf_x = [0] * (2 * Nsmooth + 1)
            psf_y = [0] * (2 * Nsmooth + 1)
            psf_z = [0] * (2 * Nsmooth + 1)
            for kk in range(2 * Nsmooth + 1):
                psf_x[kk] = 0.5 * (erfx[kk + 1] - erfx[kk])
                psf_y[kk] = 0.5 * (erfy[kk + 1] - erfy[kk])
                psf_z[kk] = 0.5 * (erfz[kk + 1] - erfz[kk])


            for kk in range(max(0, xj - Nsmooth), min(nx - 1, xj + Nsmooth)):
                for ll in range(max(0, yj - Nsmooth), min(ny - 1, yj + Nsmooth)):
                    xy_map[kk][ll] += mass[jj] * psf_x[kk - xj - Nsmooth] * psf_y[ll - yj - Nsmooth]
                for ll in range(max(0, zj - Nsmooth), min(nz - 1, zj + Nsmooth)):
                    xz_map[kk][ll] += mass[jj] * psf_x[kk - xj - Nsmooth] * psf_z[ll - zj - Nsmooth]


        xy_map *= dSxyinv
        xy_sum += xy_map
        xz_map *= dSxzinv
        xz_sum += xz_map
        xy_map = np.maximum(xy_map, fmin)
        xy_map = np.minimum(xy_map, fmax)
        xz_map = np.maximum(xz_map, fmin)
        xz_map = np.minimum(xz_map, fmax)
        # for jj in range(nx):
        #     for kk in range(ny):
        #         xy_map[jj][kk] *= dSxyinv
        #         xy_sum[jj][kk] += xy_map[jj][kk]
        #     for kk in range(nz):
        #         xz_map[jj][kk] *= dSxzinv
        #         xz_sum[jj][kk] += xz_map[jj][kk]


        # plot the data
        xz_idx = ii * nypanel + 1
        xy_idx = ii * nypanel
        img = ax[xz_idx].imshow(xz_map.T, extent = [xmin, xmax, zmin, zmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img = ax[xy_idx].imshow(xy_map.T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

        ax[xz_idx].set_ylim([zmin, zmax])
        ax[xy_idx].set_ylim([ymin, ymax])

        ax[xz_idx].spines["bottom"].set_color("white")
        ax[xy_idx].spines[   "top"].set_color("white")
        ax[xz_idx].spines[  "left"].set_color("white")
        ax[xy_idx].spines[  "left"].set_color("white")
        if ii != (nxpanel - 1):
            ax[xz_idx].spines["right"].set_color("white")
            ax[xy_idx].spines["right"].set_color("white")

        ax[xy_idx].set_xlabel(r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))

    # close the HDF5 file
    h5file.close()

    xy_sum = np.maximum(xy_sum, fmin)
    xy_sum = np.minimum(xy_sum, fmax)
    xz_sum = np.maximum(xz_sum, fmin)
    xz_sum = np.minimum(xz_sum, fmax)

    img = ax[1].imshow(xz_sum.T, extent = [xmin, xmax, zmin, zmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
    img = ax[0].imshow(xy_sum.T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

    ax[1].set_ylim([zmin, zmax])
    ax[0].set_ylim([ymin, ymax])

    ax[1].spines["bottom"].set_color("white")
    ax[0].spines[   "top"].set_color("white")
    ax[1].spines[ "right"].set_color("white")
    ax[0].spines[ "right"].set_color("white")

    ax[0].set_xlabel(r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
    ax[0].set_ylabel(r"$y$ ({:<})".format(length_unit[0].decode('UTF-8')))
    ax[1].set_ylabel(r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))



    for idx in range(nxpanel * nypanel):
        ax[idx].set_xlim([xmin, xmax])

        # ax[idx].set_xticks([-9, -6, -3, 0, 3, 6, 9])
        # ax[idx].set_yticks([-9, -6, -3, 0, 3, 6, 9])
        ax[idx].set_xticks([-10, -5, 0, 5, 10])
        ax[idx].set_yticks([-10, -5, 0, 5, 10])
        # ax[idx].grid()
        ax[idx].tick_params(axis = "both", direction = "in", color = "white", bottom = "on", top = "on", left = "on", right = "on")


        # # set caption
        # caption  = "(" + "{:^c}".format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ")"
        # caption += " " + r"${:<} = {:.1f}$ {:<}".format("z_\mathrm{d}", zd[kk], length_unit[0].decode('UTF-8'))
        # caption += ", $t =$" + r"${:.0f}$ {:<}".format(time[0] / 1000, "Gyr")
        # # ax[idx].text(xmin + 2, zmax - 7, caption, fontsize=12)
        # ax[idx].text(xmin + 0.5, zmax - 2.5, caption, fontsize=11)


    # utils.set_shared_xlabel(ax, r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
    # utils.set_shared_ylabel(ax, r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))


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
    # fig.suptitle(r"$t = {:.0f}$ {:<}".format(time[0], "Myr"))
    fig.suptitle(r"$t = {:.3f}$ {:<}".format(time[0] / 1000, "Gyr"))


    # plt.show()
    figname = "fig/" + "star" + snapshot
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


# read number of all component(s)
target = "dat/" + filename + ".profile.h5"
data_file = h5py.File(target, "r")
Nkind = data_file["/"].attrs["kinds"][0]
data_file.close()

nxpanel += Nkind - Ndark

if Nkind > Ndark:
    my_cmap = utils.generate_cmap(["darkblue", "deepskyblue", "lime", "yellow", "red", "magenta", "white"])
    cores = mp.cpu_count()
    # cores = int(np.ceil(cores / 2))
    pool = mp.Pool(cores)
    args = [(ii, Nkind) for ii in range(init, last + 1, 1)]
    pool.map(wrapper, args)
    pool.close()
