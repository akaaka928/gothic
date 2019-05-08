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


# monochrome = False
monochrome = True


nxpanel = 5
# before merger, pericentric passage, after 0, after 1, current
fileid = [0, 60, 136, 212, 288]






col_M31disk = "white"
lw_M31disk = 0.25

col_shell = "white"
# ms_shell = 3
ms_shell = 1

# col_field, col_gss, col_sc, col_sd = "white", "red", "magenta", "yellow"
col_field, col_gss, col_sc, col_sd = "white", "white", "magenta", "yellow"
lw_field, lw_gss, lw_sc, lw_sd = 0.25, 0.25, 0.25, 0.25
ms_gss, ms_sc, ms_sd = 1, 1, 1

col_wmbh = "black"
ms_wmbh = 3


pt = ["o", "s", "^", "D"]
ls = ["-", ":", "-.", "--"]
col = ["black", "red", "blue", "magenta"]


filename = "k17disk"
# Nskip = 2
# skip = [0, 2]
Nskip = 1
skip = [2]
Nmbh = 1
wmbh = [2]
fmin, fmax = 1.0e+3, 1.0e+8
xtics = [-4.0, -2.0, 0.0, 2.0, 4.0, 6.0]
ytics = [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0]
ztics = [720.0, 770.0, 820.0, 870.0]
lab = ["halo", "bulge", "MBH", "disk", "retro"]

cap = ["dark matter halo", "bulge", "prograding disk", "retrograding disk", "bulge + disks"]

col_frame = "black"
col_grid  = "white"
col_caption = "white"
if monochrome:
    col_M31disk = "black"
    col_shell = "black"
    col_field, col_gss, col_sc, col_sd = "black", "black", "black", "black"
    col_grid  = "black"
    col_caption = "black"
    # fmin, fmax = 1.0e+2, 1.0e+7


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

# read number of all component(s)
txtfile = open("doc/" + filename + ".summary.txt", "r")
unit = int(txtfile.readline())
line = txtfile.readline()
item = line.split("\t")
kind = int(item[0])
# sphe = int(item[1])
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
if monochrome:
    my_cmap = "gray_r"
    # my_cmap = "gray"


nypanel = kind
# bulge + BH, prograde disk, retrograde disk, stellar, DM
fig_xy = utils.set_figure(nxpanel, nypanel)
ax_xy = [0] * nxpanel * nypanel
utils.locate_panels(fig_xy, ax_xy, nxpanel, nypanel, True, True)
fig_xz = utils.set_figure(nxpanel, nypanel)
ax_xz = [0] * nxpanel * nypanel
utils.locate_panels(fig_xz, ax_xz, nxpanel, nypanel, True, True)

time = [0] * nxpanel

for ii in range(nxpanel):
    snapshot = "{:03d}".format(fileid[ii])
    input_file = "dat/" + filename + ".m31obs" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    time[ii] = h5file["/"].attrs["time"][0]
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

    Ndata = 0
    for jj in range(kind):
        if (Nskip == 0) or (jj not in skip):
            # read surface density maps
            folder = "field" + str(jj) + "/"
            xy_map[Ndata] = h5file[folder + "Sigma_xy"]
            xz_map[Ndata] = h5file[folder + "Sigma_zx"].value.T
            xx[Ndata] = h5file[folder +  "xi"]
            yy[Ndata] = h5file[folder + "eta"]
            zz[Ndata] = h5file[folder +   "D"]

            Ndata += 1
        if (Nmbh > 0) and (jj in wmbh):
            # read particle location
            folder = "obs" + str(jj) + "/"
            mbh_xi   = h5file[folder +   "xi"]
            mbh_eta  = h5file[folder +  "eta"]
            mbh_dist = h5file[folder + "dist"]
            mbh_obs[0] = mbh_xi[0]
            mbh_obs[1] = mbh_eta[0]
            mbh_obs[2] = mbh_dist[0]
    # close the HDF5 file
    h5file.close()

    Ndata += 1
    # all stellar components
    for jj in range(1, Ndata - 1):
        xy_map[Ndata - 1] += xy_map[jj]
        xz_map[Ndata - 1] += xz_map[jj]


    xy_map = np.maximum(xy_map, fmin)
    xy_map = np.minimum(xy_map, fmax)
    xz_map = np.maximum(xz_map, fmin)
    xz_map = np.minimum(xz_map, fmax)


    # plot the data
    xmin, xmax = xx[0][0], xx[0][nx]
    ymin, ymax = yy[0][0], yy[0][ny]
    zmin, zmax = zz[0][0], zz[0][nz]


    for jj in range(Ndata):
        img_xy = ax_xy[ii * nypanel + kind - 1 - jj].imshow(xy_map[jj].T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img_xz = ax_xz[ii * nypanel + kind - 1 - jj].imshow(xz_map[jj].T, extent = [xmin, xmax, zmin, zmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

    # incidate wandering MBH location
    if Nmbh == 1:
        # for jj in range(Ndata):
        #     ax_xy[ii * nypanel + jj].plot(mbh_obs[0], mbh_obs[1], "+", color = col_wmbh, markersize = ms_wmbh)
        #     ax_xz[ii * nypanel + jj].plot(mbh_obs[0], mbh_obs[2], "+", color = col_wmbh, markersize = ms_wmbh)
        ax_xy[ii * nypanel].plot(mbh_obs[0], mbh_obs[1], "+", color = col_wmbh, markersize = ms_wmbh)
        ax_xz[ii * nypanel].plot(mbh_obs[0], mbh_obs[2], "+", color = col_wmbh, markersize = ms_wmbh)


for idx in range(nxpanel * nypanel):
    # draw reference ellipse of the M31 disk
    ax_xy[idx].plot(disk_xi[                    ::int((Ndisk - 1) / 2)], disk_eta[                    ::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# minor axis
    ax_xy[idx].plot(disk_xi[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_eta[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# major axis
    ax_xy[idx].plot(disk_xi                                            , disk_eta                                            , "-", color = col_M31disk, linewidth = lw_M31disk)
    ax_xz[idx].plot(disk_xi[                    ::int((Ndisk - 1) / 2)], disk_D  [                    ::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# minor axis
    ax_xz[idx].plot(disk_xi[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_D  [int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# major axis
    ax_xz[idx].plot(disk_xi                                            , disk_D                                              , "-", color = col_M31disk, linewidth = lw_M31disk)

    # reference points of shells and GSS
    ax_xy[idx].plot(Eshell_xi, Eshell_eta, "o", color = col_shell, markerfacecolor = "none", markersize = ms_shell)
    ax_xy[idx].plot(Wshell_xi, Wshell_eta, "o", color = col_shell, markerfacecolor = "none", markersize = ms_shell)
    # for jj in range(Nfield):
    #     ax_xy[idx].plot(field_xi[jj], field_eta[jj], "-", color = col_field, linewidth = lw_field)

    # distance measurements to GSS
    for jj in range(Ngss):
        ax_xy[idx].plot(gss_field_xi[jj], gss_field_eta[jj], "-", color = col_gss, linewidth = lw_gss)
    ax_xz[idx].plot(gss_xi, gss_D, "s", color = col_gss, markerfacecolor = "none", markersize = ms_gss)
    ax_xz[idx].errorbar(gss_xi, gss_D, yerr = gss_Derr, ls = "none", ecolor = col_gss, elinewidth = lw_gss)


for ii in range(nxpanel * nypanel):
    ax_xy[ii].set_xlim([xmin, xmax])
    ax_xy[ii].set_ylim([ymin, ymax])
    ax_xy[ii].set_xticks(xtics)
    ax_xy[ii].set_yticks(ytics)
    ax_xy[ii].tick_params(axis = "both", direction = "in", color = col_grid, bottom = True, top = True, left = True, right = True)
    ax_xy[ii].spines["bottom"].set_color(col_grid)
    ax_xy[ii].spines[   "top"].set_color(col_grid)
    ax_xy[ii].spines[  "left"].set_color(col_grid)
    ax_xy[ii].spines[ "right"].set_color(col_grid)

    ax_xz[ii].set_xlim([xmin, xmax])
    ax_xz[ii].set_ylim([zmin, zmax])
    ax_xz[ii].set_xticks(xtics)
    ax_xz[ii].set_yticks(ztics)
    ax_xz[ii].tick_params(axis = "both", direction = "in", color = col_grid, bottom = True, top = True, left = True, right = True)
    ax_xz[ii].spines["bottom"].set_color(col_grid)
    ax_xz[ii].spines[   "top"].set_color(col_grid)
    ax_xz[ii].spines[  "left"].set_color(col_grid)
    ax_xz[ii].spines[ "right"].set_color(col_grid)

# set caption
t0 = time[nxpanel - 1]
time -= t0
for ii in range(nxpanel):
    for jj in range(nypanel):
        caption  = "(" + "{:^c}".format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ")"
        if ii == 0:
            caption += cap[nypanel - 1 - jj]
        idx = ii * nypanel + jj
        ax_xy[idx].text(xmin + 0.03 * (xmax - xmin), ymax - 0.12 * (ymax - ymin), caption, color = col_caption, fontsize = 11)
        ax_xz[idx].text(xmin + 0.03 * (xmax - xmin), zmax - 0.12 * (zmax - zmin), caption, color = col_caption, fontsize = 11)
        if jj == 0:
            ax_xy[idx].text(xmax - 0.05 * (xmax - xmin), ymin + 0.08 * (ymax - ymin), r"$t = {:.1f}$ {:<}".format(time[ii], "Myr"), color = col_caption, fontsize = 11, ha = "right")
            ax_xz[idx].text(xmax - 0.05 * (xmax - xmin), zmin + 0.08 * (zmax - zmin), r"$t = {:.1f}$ {:<}".format(time[ii], "Myr"), color = col_caption, fontsize = 11, ha = "right")



for ii in range(nxpanel):
    ax_xy[ii * nypanel].set_xlabel(r"$\xi$ ({:<})".format("degree"))
    ax_xy[ii * nypanel + nypanel - 1].spines[   "top"].set_color(col_frame)
    ax_xy[ii * nypanel].spines["bottom"].set_color(col_frame)
    ax_xz[ii * nypanel].set_xlabel(r"$\xi$ ({:<})".format("degree"))
    ax_xz[ii * nypanel + nypanel - 1].spines[   "top"].set_color(col_frame)
    ax_xz[ii * nypanel].spines["bottom"].set_color(col_frame)

for ii in range(nypanel):
    ax_xy[ii].set_ylabel(r"$\eta$ ({:<})".format("degree"))
    ax_xy[ii].spines["left"].set_color(col_frame)
    ax_xy[(nxpanel - 1) * nypanel + ii].spines["right"].set_color(col_frame)
    ax_xz[ii].set_ylabel(r"$D$ ({:<})".format("kpc"))
    ax_xz[ii].spines["left"].set_color(col_frame)
    ax_xz[(nxpanel - 1) * nypanel + ii].spines["right"].set_color(col_frame)


# add colorbar
x0, x1 = 1, 0
y0, y1 = 1, 0
for at in ax_xy:
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
colorbar_ax = fig_xy.add_axes([x1, y0, 0.1 / nxpanel, y1 - y0])
cbar = fig_xy.colorbar(img_xy, cax = colorbar_ax, label = r"$\Sigma$ ($M_\odot$ deg$^{-2}$)")
cbar.solids.set_edgecolor("face")

x0, x1 = 1, 0
y0, y1 = 1, 0
for at in ax_xz:
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
colorbar_ax = fig_xz.add_axes([x1, y0, 0.1 / nxpanel, y1 - y0])
cbar = fig_xz.colorbar(img_xz, cax = colorbar_ax, label = r"$\Sigma$ ($M_\odot$ deg$^{-2}$)")
cbar.solids.set_edgecolor("face")


# # add current time
# fig.suptitle(r"$t = {:.2f}$ {:<}".format(time, "Myr"))
# # fig.suptitle(r"$t = {:.3f}$ {:<}".format(time / 1000, "Gyr"))


# save figures
# figname = "fig/" + filename + "_map" + snapshot
figname = filename + "_evol"
if monochrome:
    figname += "_mono"
fig_xy.savefig(figname + "_xy" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
fig_xz.savefig(figname + "_xz" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
if outputPDF:
    fig_xy.savefig(figname + "_xy" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    fig_xz.savefig(figname + "_xz" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.close("all")
