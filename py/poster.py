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


outputPDF = True
# outputPDF = False


col_M31disk = "white"
# lw_M31disk = 0.5
lw_M31disk = 1

col_shell = "white"
# ms_shell = 3
ms_shell = 8

col_field, col_gss, col_sc, col_sd = "white", "red", "magenta", "yellow"
# lw_field, lw_gss, lw_sc, lw_sd = 0.5, 0.5, 0.5, 0.5
lw_field, lw_gss, lw_sc, lw_sd = 2, 2, 2, 2
ms_gss, ms_sc, ms_sd = 3, 3, 3

col_wmbh = "black"
# ms_wmbh = 3
ms_wmbh = 16


pt = ["o", "s", "^", "D"]
ls = ["-", ":", "-.", "--"]
col = ["black", "red", "blue", "magenta"]


filename = "k17disk"
Nskip = 2
skip = [0, 2]
Nmbh = 1
wmbh = [2]
init = 0
last = 399
# init = 270
# last = 270
fmin, fmax = 1.0e+4, 1.0e+9
xtics = [-4.0, -2.0, 0.0, 2.0, 4.0, 6.0]
ytics = [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0]
ztics = [700.0, 750.0, 800.0, 850.0, 900.0]
lab = ["halo", "bulge", "MBH", "disk", "retro"]


def draw_figure(fileid, kind, Ndisk, disk_xi, disk_eta, disk_D, Neast, Eshell_xi, Eshell_eta, Nwest, Wshell_xi, Wshell_eta, Nfield, field_xi, field_eta, Ngss, gss_xi, gss_eta, gss_D, gss_Derr, gss_field_xi, gss_field_eta, NstreamC, streamC_xi, streamC_eta, streamC_D, streamC_Derr, streamC_field_xi, streamC_field_eta, NstreamD, streamD_xi, streamD_eta, streamD_D, streamD_Derr, streamD_field_xi, streamD_field_eta):
    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".m31obs" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    time = h5file["/"].attrs["time"][0]

    # kind = h5file["/"].attrs["kinds"][0]
    nx = h5file["/"].attrs["nx"][0]
    ny = h5file["/"].attrs["ny"][0]
    nz = h5file["/"].attrs["nz"][0]

    # memory allocation for surface density maps
    xy_map = np.zeros((kind, nx, ny))
    xz_map = np.zeros((kind, nx, nz))
    xx     = np.zeros((kind, nx + 1))
    yy     = np.zeros((kind, ny + 1))
    zz     = np.zeros((kind, nz + 1))

    mbh_obs = [0] * 3

    nxpanel = 0
    nypanel = 2
    Ndata = 0
    for ii in range(kind):
        if (Nskip == 0) or (ii not in skip):
            # read surface density maps
            folder = "field" + str(ii) + "/"
            xy_map[nxpanel] = h5file[folder + "Sigma_xy"]
            xz_map[nxpanel] = h5file[folder + "Sigma_zx"][()].T
            xx[nxpanel] = h5file[folder +  "xi"]
            yy[nxpanel] = h5file[folder + "eta"]
            zz[nxpanel] = h5file[folder +   "D"]

            nxpanel += 1
            Ndata += 1
        if (Nmbh > 0) and (ii in wmbh):
            # read particle position and mass
            folder = "obs" + str(ii) + "/"
            mbh_xi   = h5file[folder +   "xi"]
            mbh_eta  = h5file[folder +  "eta"]
            mbh_dist = h5file[folder + "dist"]
            mbh_obs[0] = mbh_xi[0]
            mbh_obs[1] = mbh_eta[0]
            mbh_obs[2] = mbh_dist[0]
    # close the HDF5 file
    h5file.close()

    summary = False
    if nxpanel > 1:
        for ii in range(nxpanel):
            xy_map[nxpanel] += xy_map[ii]
            xz_map[nxpanel] += xz_map[ii]

        nxpanel += 1
        summary = True


    xy_map = np.maximum(xy_map, fmin)
    xy_map = np.minimum(xy_map, fmax)
    xz_map = np.maximum(xz_map, fmin)
    xz_map = np.minimum(xz_map, fmax)


    # plot the data
    nxpanel_plt = 1
    nypanel_plt = 1
    fig = utils.set_figure(nxpanel_plt, nypanel_plt)
    ax = [0] * nxpanel_plt * nypanel_plt
    utils.locate_panels(fig, ax, nxpanel_plt, nypanel_plt, True, True)

    xmin, xmax = xx[0][0], xx[0][nx]
    ymin, ymax = yy[0][0], yy[0][ny]
    zmin, zmax = zz[0][0], zz[0][nz]

    head = 0
    if summary:
        # img = ax[1].imshow(xz_map[nxpanel - 1].T, extent = [xmin, xmax, zmin, zmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img = ax[0].imshow(xy_map[nxpanel - 1].T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        head = 1

    # for ii in range(Ndata):
    #     img = ax[(head + ii) * nypanel + 1].imshow(xz_map[ii].T, extent = [xmin, xmax, zmin, zmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
    #     img = ax[(head + ii) * nypanel    ].imshow(xy_map[ii].T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")



    for ii in range(nxpanel_plt):
        idx = ii * nypanel_plt
        # reference ellipse of the M31 disk
        ax[idx    ].plot(disk_xi[                    ::int((Ndisk - 1) / 2)], disk_eta[                    ::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# minor axis
        ax[idx    ].plot(disk_xi[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_eta[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# major axis
        ax[idx    ].plot(disk_xi                                            , disk_eta                                            , "-", color = col_M31disk, linewidth = lw_M31disk)
        # ax[idx + 1].plot(disk_xi[                    ::int((Ndisk - 1) / 2)], disk_D  [                    ::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# minor axis
        # ax[idx + 1].plot(disk_xi[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_D  [int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# major axis
        # ax[idx + 1].plot(disk_xi                                            , disk_D                                              , "-", color = col_M31disk, linewidth = lw_M31disk)


        # reference points of shells and GSS
        ax[idx].plot(Eshell_xi, Eshell_eta, "o", color = col_shell, markerfacecolor = "none", markersize = ms_shell)
        ax[idx].plot(Wshell_xi, Wshell_eta, "o", color = col_shell, markerfacecolor = "none", markersize = ms_shell)
        # for jj in range(Nfield):
        #     ax[idx].plot(field_xi[jj], field_eta[jj], "-", color = col_field, linewidth = lw_field)

        # distance measurements to GSS
        for jj in range(Ngss):
            ax[idx].plot(gss_field_xi[jj], gss_field_eta[jj], "-", color = col_gss, linewidth = lw_gss)
        # ax[idx + 1].plot(gss_xi, gss_D, "s", color = col_gss, markerfacecolor = "none", markersize = ms_gss)
        # ax[idx + 1].errorbar(gss_xi, gss_D, yerr = gss_Derr, ls = "none", ecolor = col_gss, elinewidth = lw_gss)

        # # distance measurements to Stream C
        # for jj in range(NstreamC):
        #     ax[idx].plot(streamC_field_xi[jj], streamC_field_eta[jj], "-", color = col_sc, linewidth = lw_sc)
        # # ax[idx + 1].plot(streamC_xi, streamC_D, "s", color = col_sc, markerfacecolor = "none", markersize = ms_sc)
        # # ax[idx + 1].errorbar(streamC_xi, streamC_D, yerr = streamC_Derr, ls = "none", ecolor = col_sc, elinewidth = lw_sc)

        # distance measurements to Stream D
        for jj in range(NstreamD):
            ax[idx].plot(streamD_field_xi[jj], streamD_field_eta[jj], "-", color = col_sd, linewidth = lw_sd)
        # ax[idx + 1].plot(streamD_xi, streamD_D, "s", color = col_sd, markerfacecolor = "none", markersize = ms_sd)
        # ax[idx + 1].errorbar(streamD_xi, streamD_D, yerr = streamD_Derr, ls = "none", ecolor = col_sd, elinewidth = lw_sd)


    # incidate wandering MBH location
    if Nmbh == 1:
        ax[0].plot(mbh_obs[0], mbh_obs[1], "+", color = col_wmbh, markersize = ms_wmbh)
        # ax[1].plot(mbh_obs[0], mbh_obs[2], "+", color = col_wmbh, markersize = ms_wmbh)


    for ii in range(nxpanel_plt * nypanel_plt):
        ax[ii].set_xlim([xmin, xmax])
        ax[ii].set_xticks(xtics)
        # ax[ii].set_yticks(ytics)
        ax[ii].tick_params(axis = "both", direction = "in", color = "white", bottom = True, top = True, left = True, right = True)
        ax[ii].spines["bottom"].set_color("white")
        ax[ii].spines[   "top"].set_color("white")
        ax[ii].spines[  "left"].set_color("white")
        ax[ii].spines[ "right"].set_color("white")

    for ii in range(nxpanel_plt):
        # ax[ii * nypanel_plt + 1].set_yticks(ztics)
        ax[ii * nypanel_plt    ].set_yticks(ytics)
        # ax[ii * nypanel_plt + 1].set_ylim([zmin, zmax])
        ax[ii * nypanel_plt    ].set_ylim([ymin, ymax])
        ax[ii * nypanel_plt    ].set_xlabel(r"$\xi$ ({:<})".format("degree"))
        # ax[ii * nypanel_plt + 1].spines[   "top"].set_color("black")
        ax[ii * nypanel_plt    ].spines["bottom"].set_color("black")


    ax[0].set_ylabel(r"$\eta$ ({:<})".format("degree"))
    # ax[1].set_ylabel(r"$D$ ({:<})".format("kpc"))
    ax[0].spines["left"].set_color("black")
    # ax[1].spines["left"].set_color("black")
    # ax[(nxpanel_plt - 1) * nypanel_plt + 1].spines["right"].set_color("black")
    ax[(nxpanel_plt - 1) * nypanel_plt    ].spines["right"].set_color("black")


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
    # colorbar_ax = fig.add_axes([x1, y0, 0.1 / nxpanel_plt, y1 - y0])
    colorbar_ax = fig.add_axes([x1, y0, 0.05 / nxpanel_plt, y1 - y0])
    cbar = fig.colorbar(img, cax = colorbar_ax, label = r"$\Sigma$ ($M_\odot$ deg$^{-2}$)")
    cbar.solids.set_edgecolor("face")

    # # add current time
    # fig.suptitle(r"$t = {:.2f}$ {:<}".format(time, "Myr"))
    # # fig.suptitle(r"$t = {:.3f}$ {:<}".format(time / 1000, "Gyr"))


    # save figures
    figname = "fig/" + filename + "_poster" + snapshot
    fig.savefig(figname + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if outputPDF:
        fig.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.close("all")



def wrapper(argv):
    return draw_figure(*argv)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
plt.rcParams['font.size'] = 28
# plt.rcParams['font.size'] = 14

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# read number of all component(s)
txtfile = open("doc/" + filename + ".summary.txt", "r")
unit = int(txtfile.readline())
line = txtfile.readline()
item = line.split("\t")
Nkind = int(item[0])
# Nsphe = int(item[1])
txtfile.close()

# set reference ellipse of M31 disk
Ndisk, disk_xi, disk_eta, disk_D = m31.disk_ellipse()

# set reference points of shells and GSS
Neast, Eshell_xi, Eshell_eta, Nwest, Wshell_xi, Wshell_eta = m31.Fardal2007_shell()

# set reference points of distance measurements by Conn et al. (2016)
Nfield, field_xi, field_eta = m31.GSS_obs_field()
Ngss, gss_xi, gss_eta, gss_D, gss_Dep, gss_Dem, gss_field_xi, gss_field_eta = m31.GSS_distance()
NstreamC, streamC_xi, streamC_eta, streamC_D, streamC_Dep, streamC_Dem, streamC_field_xi, streamC_field_eta = m31.StreamC_distance()
NstreamD, streamD_xi, streamD_eta, streamD_D, streamD_Dep, streamD_Dem, streamD_field_xi, streamD_field_eta = m31.StreamD_distance()
gss_Derr = np.zeros((2, Ngss))
for ii in range(Ngss):
    gss_Derr[0][ii] = gss_Dem[ii]
    gss_Derr[1][ii] = gss_Dep[ii]
streamC_Derr = np.zeros((2, NstreamC))
for ii in range(NstreamC):
    streamC_Derr[0][ii] = streamC_Dem[ii]
    streamC_Derr[1][ii] = streamC_Dep[ii]
streamD_Derr = np.zeros((2, NstreamD))
for ii in range(NstreamD):
    streamD_Derr[0][ii] = streamD_Dem[ii]
    streamD_Derr[1][ii] = streamD_Dep[ii]


my_cmap = utils.generate_cmap(["darkblue", "deepskyblue", "lime", "yellow", "red", "magenta", "white"])
# cores = mp.cpu_count()
cores = int(np.ceil(mp.cpu_count() / 2))
pool = mp.Pool(cores)
args = [(ii, Nkind, Ndisk, disk_xi, disk_eta, disk_D, Neast, Eshell_xi, Eshell_eta, Nwest, Wshell_xi, Wshell_eta, Nfield, field_xi, field_eta, Ngss, gss_xi, gss_eta, gss_D, gss_Derr, gss_field_xi, gss_field_eta, NstreamC, streamC_xi, streamC_eta, streamC_D, streamC_Derr, streamC_field_xi, streamC_field_eta, NstreamD, streamD_xi, streamD_eta, streamD_D, streamD_Derr, streamD_field_xi, streamD_field_eta) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
