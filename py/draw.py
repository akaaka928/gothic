import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm # for logarithmic plot in imshow
import matplotlib.ticker as ticker

import os.path as path


from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# import multiprocessing as mp

import utils as utils


outputPDF = False

fs_base = 24
tl_base = 6.0
tw_base = 1.0


# filename = "mw"
# Nskip = 1
# skip = [0]
# init = 0
# last = 127
# # last = 0
# fmin, fmax = 1.0e-3, 1.0e+1
# fvmin, fvmax = 1.0e+4, 1.0e+8
# xtics = [-10, -5, 0, 5, 10]
# ytics = [-10, -5, 0, 5, 10]
# lab = ["DM halo", "stellar halo", "bulge", "central BH", "thick disk", "thin disk"]

filename = "m31"
Nskip = 1
skip = [0]
init = 0
last = 47
fmin, fmax = 1.0e-3, 1.0e+1
fvmin, fvmax = 1.0e+4, 1.0e+8
xtics = [-10, -5, 0, 5, 10]
ytics = [-10, -5, 0, 5, 10]
lab = ["DM halo", "stellar halo", "bulge", "disk"]

# filename = "k17disk"
# Nskip = 2
# skip = [0, 2]
# init = 0
# last = 399
# fmin, fmax = 1.0e-7, 1.0e-1
# fvmin, fvmax = 1.0e+0, 1.0e+6
# lab = ["halo", "bulge", "MBH", "disk"]

# # filename = "hd"
# filename = "disk"
# Nskip = 0
# init = 0
# last = 13
# fmin, fmax = 1.0e-3, 1.0e+0
# fvmin, fvmax = 1.0e-3, 1.0e+0
# lab = ["halo", "disk"]

# # filename = "halocusp_run"
# # filename = "halocore1_run"
# filename = "halocore2_run"
# # filename = "halocore3_run"
# Nskip = 1
# skip = [0]
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-6, 1.0e-1
# fvmin, fvmax = 1.0e-6, 1.0e-1
# lab = ["halo", "core", "GC"]

# filename = "hernquist"
# Nskip = 0
# init = 0
# last = 47
# # last = 0
# fmin, fmax = 1.0e-3, 1.0e+1
# fvmin, fvmax = 1.0e-3, 1.0e+1
# lab = ["halo"]

# filename = "etg"
# Nskip = 1
# skip = [0]
# # Nskip = 0
# init = 0
# last = 47
# last = 15
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["DM halo", "bulge", "stellar halo", "central BH"]

# filename = "satellite"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["bulge"]

# filename = "w3iso"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w3ra2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w3ra1"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w3ra1_2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w3ra1_4"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w3ra1_8"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w5iso"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w5ra2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w5ra1"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w5ra1_2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w5ra1_4"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w5ra1_8"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w7iso"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w7ra2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w7ra1"
# Nskip = 0
# init = 0
# last = 140
# last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w7ra1_2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w7ra1_4"
# Nskip = 0
# init = 0
# last = 140
# last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "w7ra1_8"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+3, 1.0e+6
# lab = ["GC"]

# filename = "m12iso"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 3.1e+4, 1.0e+8
# lab = ["halo"]

# filename = "m12ra2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 3.1e+4, 1.0e+8
# lab = ["halo"]

# filename = "m12ra1"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 3.1e+4, 1.0e+8
# lab = ["halo"]

# filename = "m12ra1_2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 3.1e+4, 1.0e+8
# lab = ["halo"]

# filename = "m12ra1_4"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 3.1e+4, 1.0e+8
# lab = ["halo"]

# filename = "m12ra1_8"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 3.1e+4, 1.0e+8
# lab = ["halo"]

# filename = "m09iso"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 3.1e-2
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "m09ra2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 3.1e-2
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "m09ra1"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 3.1e-2
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "m09ra1_2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 3.1e-2
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "m09ra1_4"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 3.1e-2
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "m09ra1_8"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 3.1e-2
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "a00iso"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "a00ra2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "a00ra1"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "a00ra1_2"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]

# filename = "a00ra1_4"
# Nskip = 0
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-4, 1.0e-1
# fvmin, fvmax = 1.0e+4, 3.1e+7
# lab = ["halo"]


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", ":", "-.", "--"]
Ncol = 8
col = ["black", "red", "blue", "magenta", "green", "brown", "cyan", "yellow"]

pt_mbh = "+"
col_mbh = "black"
ms_mbh_base = 8


def draw_figure(fileid, kind, sphe):
    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".plt" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    length_unit = h5file["/"].attrs["length_astro_unit_name"][0]
    time_unit = h5file["/"].attrs["time_astro_unit_name"][0]
    mass_unit = h5file["/"].attrs["mass_astro_unit_name"][0]
    velocity_unit = h5file["/"].attrs["velocity_astro_unit_name"][0]
    col_density_unit = h5file["/"].attrs["col_density_astro_unit_name"][0]
    time = h5file["/"].attrs["time"][0]

    nx = h5file["/"].attrs["nx"][0]
    ny = h5file["/"].attrs["ny"][0]
    nz = h5file["/"].attrs["nz"][0]
    nv = h5file["/"].attrs["nz"][0]

    # memory allocation for surface density maps
    xy_map = np.zeros((kind + 1, nx, ny))
    zy_map = np.zeros((kind + 1, nz, ny))
    xx     = np.zeros((kind + 1, nx + 1))
    yy     = np.zeros((kind + 1, ny + 1))
    zz     = np.zeros((kind + 1, nz + 1))
    xv_map = np.zeros((kind + 1, nx, nv))
    yv_map = np.zeros((kind + 1, ny, nv))
    zv_map = np.zeros((kind + 1, nz, nv))
    vv     = np.zeros((kind + 1, nv + 1))

    # memory allocation for MBH location
    mbh_x = [0] * sphe
    mbh_y = [0] * sphe
    mbh_z = [0] * sphe
    mbh_vx = [0] * sphe
    mbh_vy = [0] * sphe
    mbh_vz = [0] * sphe
    Nmbh = 0


    nxpanel = 2
    nypanel = 0
    nvpanel = 3
    Ndata = 0
    for ii in range(kind):
        num = h5file["attr" + str(ii)].attrs["number"][0]
        if (num > 1):
            if (Nskip == 0) or (ii not in skip):
                # read surface density maps
                folder = "field" + str(ii) + "/"
                xy_map[nypanel] = h5file[folder + "Sigma_xy"]
                zy_map[nypanel] = h5file[folder + "Sigma_yz"][()].T
                xv_map[nypanel] = h5file[folder + "f_xv"]
                yv_map[nypanel] = h5file[folder + "f_yv"]
                zv_map[nypanel] = h5file[folder + "f_zv"]
                xx[nypanel] = h5file[folder + "x"]
                yy[nypanel] = h5file[folder + "y"]
                zz[nypanel] = h5file[folder + "z"]
                vv[nypanel] = h5file[folder + "v"]

                nypanel += 1
                Ndata += 1

        else:
            # read MBH location
            # tmpfile = h5py.File("dat/" + filename + ".split" + snapshot + ".h5", "r")
            # position = tmpfile["data" + str(ii) + "/position"]
            position = h5file["/attr" + str(ii)].attrs["com"]
            velocity = h5file["/attr" + str(ii)].attrs["vel"]
            mbh_x[Nmbh] = position[0]
            mbh_y[Nmbh] = position[1]
            mbh_z[Nmbh] = position[2]
            mbh_vx[Nmbh] = velocity[0]
            mbh_vy[Nmbh] = velocity[1]
            mbh_vz[Nmbh] = velocity[2]
            Nmbh += 1


    # close the HDF5 file
    h5file.close()

    summary = False
    if nypanel > 1:
        for ii in range(nypanel):
            xy_map[nypanel] += xy_map[ii]
            zy_map[nypanel] += zy_map[ii]
            xv_map[nypanel] += xv_map[ii]
            yv_map[nypanel] += yv_map[ii]
            zv_map[nypanel] += zv_map[ii]

        nypanel += 1
        summary = True


    xy_map = np.maximum(xy_map, fmin)
    xy_map = np.minimum(xy_map, fmax)
    zy_map = np.maximum(zy_map, fmin)
    zy_map = np.minimum(zy_map, fmax)
    xv_map = np.maximum(xv_map, fvmin)
    xv_map = np.minimum(xv_map, fvmax)
    yv_map = np.maximum(yv_map, fvmin)
    yv_map = np.minimum(yv_map, fvmax)
    zv_map = np.maximum(zv_map, fvmin)
    zv_map = np.minimum(zv_map, fvmax)


    # plot the data
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)
    fig_v = utils.set_figure(nvpanel, nypanel)
    ax_v = [0] * nvpanel * nypanel
    utils.locate_panels(fig_v, ax_v, nvpanel, nypanel, True, True)

    xmin, xmax = xx[0][0], xx[0][nx]
    ymin, ymax = yy[0][0], yy[0][ny]
    zmin, zmax = zz[0][0], zz[0][nz]
    vmin, vmax = vv[0][0], vv[0][nv]

    # adjust marker size and etcetra
    fs = fs_base / np.sqrt(max([nxpanel, nypanel]))
    tl = tl_base / np.sqrt(max([nxpanel, nypanel]))
    tw = tw_base / np.sqrt(max([nxpanel, nypanel]))
    ms_mbh = ms_mbh_base / max([nxpanel, nypanel])
    fs_v = fs_base / np.sqrt(max([nvpanel, nypanel]))
    tl_v = tl_base / np.sqrt(max([nvpanel, nypanel]))
    tw_v = tw_base / np.sqrt(max([nvpanel, nypanel]))
    ms_mbh_v = ms_mbh_base / max([nvpanel, nypanel])

    head = 0
    if summary:
        img = ax[    nypanel - 1].imshow(xy_map[nypanel - 1].T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img = ax[2 * nypanel - 1].imshow(zy_map[nypanel - 1].T, extent = [zmin, zmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img_v = ax_v[    nypanel - 1].imshow(xv_map[nypanel - 1].T, extent = [xmin, xmax, vmin, vmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fvmin, vmax = fvmax), cmap = "hot", aspect = "auto")
        img_v = ax_v[2 * nypanel - 1].imshow(yv_map[nypanel - 1].T, extent = [ymin, ymax, vmin, vmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fvmin, vmax = fvmax), cmap = "hot", aspect = "auto")
        img_v = ax_v[3 * nypanel - 1].imshow(zv_map[nypanel - 1].T, extent = [zmin, zmax, vmin, vmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fvmin, vmax = fvmax), cmap = "hot", aspect = "auto")
        if Nmbh > 0:
            ax[    nypanel - 1].plot(mbh_x[:Nmbh], mbh_y[:Nmbh], pt_mbh, color = col_mbh, markersize = ms_mbh)
            ax[2 * nypanel - 1].plot(mbh_z[:Nmbh], mbh_y[:Nmbh], pt_mbh, color = col_mbh, markersize = ms_mbh)
            ax_v[    nypanel - 1].plot(mbh_x[:Nmbh], mbh_vx[:Nmbh], pt_mbh, color = col_mbh, markersize = ms_mbh_v)
            ax_v[2 * nypanel - 1].plot(mbh_y[:Nmbh], mbh_vy[:Nmbh], pt_mbh, color = col_mbh, markersize = ms_mbh_v)
            ax_v[3 * nypanel - 1].plot(mbh_z[:Nmbh], mbh_vz[:Nmbh], pt_mbh, color = col_mbh, markersize = ms_mbh_v)
        head = 1

    for ii in range(Ndata):
        img = ax[    nypanel - 1 - (head + ii)].imshow(xy_map[ii].T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img = ax[2 * nypanel - 1 - (head + ii)].imshow(zy_map[ii].T, extent = [zmin, zmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img_v = ax_v[    nypanel - 1 - (head + ii)].imshow(xv_map[ii].T, extent = [xmin, xmax, vmin, vmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fvmin, vmax = fvmax), cmap = "hot", aspect = "auto")
        img_v = ax_v[2 * nypanel - 1 - (head + ii)].imshow(yv_map[ii].T, extent = [ymin, ymax, vmin, vmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fvmin, vmax = fvmax), cmap = "hot", aspect = "auto")
        img_v = ax_v[3 * nypanel - 1 - (head + ii)].imshow(zv_map[ii].T, extent = [zmin, zmax, vmin, vmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fvmin, vmax = fvmax), cmap = "hot", aspect = "auto")

    for ii in range(nxpanel * nypanel):
        ax[ii].set_ylim([ymin, ymax])
        ax[ii].tick_params(axis = "both", direction = "in", color = "white", bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
        ax[ii].spines["bottom"].set_color("white")
        ax[ii].spines[   "top"].set_color("white")
        ax[ii].spines[  "left"].set_color("white")
        ax[ii].spines[ "right"].set_color("white")

    for jj in range(nypanel):
        ax[          jj].set_xlim([xmin, xmax])
        ax[nypanel + jj].set_xlim([zmin, zmax])
        ax[          jj].set_xlabel(r"$x$" + r" ({:<})".format(length_unit.decode("UTF-8")), fontsize = fs)
        ax[nypanel + jj].set_xlabel(r"$z$" + r" ({:<})".format(length_unit.decode("UTF-8")), fontsize = fs)
        ax[          jj].spines[ "left"].set_color("black")
        ax[nypanel + jj].spines["right"].set_color("black")
        ax[jj].set_ylabel(r"$y$" + r" ({:<})".format(length_unit.decode("UTF-8")), fontsize = fs)

    for ii in range(nxpanel):
        ax[ ii      * nypanel    ].spines["bottom"].set_color("black")
        ax[(ii + 1) * nypanel - 1].spines[   "top"].set_color("black")


    ax[0          ].set_xlabel(r"$x$" + r"~({:<})".format(length_unit.decode("UTF-8")), fontsize = fs)
    ax[    nypanel].set_xlabel(r"$z$" + r"~({:<})".format(length_unit.decode("UTF-8")), fontsize = fs)

    for ii in range(nvpanel * nypanel):
        ax_v[ii].set_ylim([vmin, vmax])
        ax_v[ii].tick_params(axis = "both", direction = "in", color = "white", bottom = True, top = True, left = True, right = True, labelsize = fs_v, length = tl_v, width = tw_v)
        ax_v[ii].spines["bottom"].set_color("white")
        ax_v[ii].spines[   "top"].set_color("white")
        ax_v[ii].spines[  "left"].set_color("white")
        ax_v[ii].spines[ "right"].set_color("white")

    for jj in range(nypanel):
        # ax_v[              jj].set_ylabel(r"$v$" + r" ({:<})".format(velocity_unit.decode("UTF-8")), fontsize = fs_v)
        ax_v[              jj].set_ylabel(r"$v$" + r"~(\si{km.s^{-1}})", fontsize = fs_v)
        ax_v[              jj].spines["left"].set_color("black")
        ax_v[              jj].set_xlim([xmin, xmax])
        ax_v[    nypanel + jj].set_xlim([ymin, ymax])
        ax_v[2 * nypanel + jj].set_xlim([zmin, zmax])
        ax_v[2 * nypanel + jj].spines["right"].set_color("black")

    for ii in range(nvpanel):
        ax_v[ ii      * nypanel    ].spines["bottom"].set_color("black")
        ax_v[(ii + 1) * nypanel - 1].spines["top"].set_color("black")

    ax_v[0          ].set_xlabel(r"$x$" + r" ({:<})".format(length_unit.decode("UTF-8")), fontsize = fs_v)
    ax_v[    nypanel].set_xlabel(r"$y$" + r" ({:<})".format(length_unit.decode("UTF-8")), fontsize = fs_v)
    ax_v[2 * nypanel].set_xlabel(r"$z$" + r" ({:<})".format(length_unit.decode("UTF-8")), fontsize = fs_v)



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
    colorbar_ax = fig.add_axes([x1, y0, 0.05 / nxpanel, y1 - y0])
    # cbar = fig.colorbar(img, cax = colorbar_ax, format=ticker.FuncFormatter(utils.scientific), label = r"$\Sigma$" + r" ({:<})".format(col_density_unit.decode("UTF-8")))
    cbar = fig.colorbar(img, cax = colorbar_ax, format=ticker.FuncFormatter(utils.scientific), label = r"$\Sigma$" + r"~(\si{g.cm^{-2}})")
    cbar.solids.set_edgecolor("face")
    # cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=fs)
    # cbar.ax.set_ylabel(cbar.ax.get_ylabel(), fontsize=fs)

    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax_v:
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
    colorbar_ax = fig_v.add_axes([x1, y0, 0.05 / nvpanel, y1 - y0])
    # cbar = fig_v.colorbar(img_v, cax = colorbar_ax, format=ticker.FuncFormatter(utils.scientific), label = r"$f$ (" + r"{:<}".format(mass_unit.decode("UTF-8")) + r" {:<}".format(length_unit.decode("UTF-8")) + r"$^{-1}$" + r" {:<}".format(velocity_unit.decode("UTF-8")) + r"$^{-1}$" + r")")
    cbar = fig_v.colorbar(img_v, cax = colorbar_ax, format=ticker.FuncFormatter(utils.scientific), label = r"$f$" + r"~(\si{M_\odot.kpc^{-1}.km^{-1}.s})")
    cbar.solids.set_edgecolor("face")
    # cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=fs_v)
    # cbar.ax.set_ylabel(cbar.ax.get_ylabel(), fontsize=fs_v)



    # add current time
    # fig.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")), fontsize = fs)
    # fig_v.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")), fontsize = fs_v)
    fig.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})", fontsize = fs)
    fig_v.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})", fontsize = fs)


    # save figures
    figname = "fig/" + filename + "_map" + snapshot
    fig.savefig(figname + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if outputPDF:
        fig.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    figname = "fig/" + filename + "_vmap" + snapshot
    fig_v.savefig(figname + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if outputPDF:
        fig_v.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.close("all")



# def wrapper(argv):
#     return draw_figure(*argv)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath,siunitx}"

# set font size
plt.rcParams['font.size'] = 18
# plt.rcParams['font.size'] = 16

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# read number of all component(s)
txtfile = open("doc/" + filename + ".summary.txt", "r")
unit = int(txtfile.readline())
line = txtfile.readline()
item = line.split("\t")
Nkind = int(item[0])
Nsphe = int(item[1])
txtfile.close()


my_cmap = utils.generate_cmap(["darkblue", "deepskyblue", "lime", "yellow", "red", "magenta", "white"])
# cores = int(np.ceil(mp.cpu_count() / 2))
# pool = mp.Pool(cores)
# args = [(ii, Nkind, Nsphe) for ii in range(init, last + 1, 1)]
# pool.map(wrapper, args)
# pool.close()
# def draw_figure(fileid, kind, sphe):
for fileid in range(init + rank, last + 1, size):
    draw_figure(fileid, Nkind, Nsphe)
