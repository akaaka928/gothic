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


monochrome = False
# monochrome = True

fs_base = 16
tl_base = 6.0
tw_base = 1.0

col_ref = "gray"
lw_ref_base = 1.0

col_M31disk = "white"
lw_M31disk_base = 1.0

col_wmbh = "black"
ms_wmbh_base = 3


pt = ["o", "s", "^", "D"]
ls = ["-", ":", "-.", "--"]
col = ["black", "red", "blue", "magenta"]

# col_hsc, col_nws, col_gcs = "white", "orange", "cyan"
col_hsc, col_nws, col_gcs = "white", "brown", "cyan"
lw_hsc_base, lw_nws_base, lw_gcs_base = 1.5, 1.5, 1.5
ms_hsc_base, ms_nws_base, ms_gcs_base = 1, 4, 4

col_frame = "black"
col_grid  = "white"
col_caption = "white"
if monochrome:
    col_M31disk = "black"
    col_grid  = "black"
    col_caption = "black"
    col_hsc, col_nws, col_gcs = "black", "black", "black"



filename = "k18nws"
Nskip = 0
Nmbh = 0
init = 0
last = 560
init = 325
last = 325
Sigma_min, Sigma_max = 3.1e+4, 1.0e+8
depth_min, depth_max = 3.1e+3, 3.1e+6
phase_min, phase_max = 1.0e+4, 3.1e+6
xtics = [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0]
ytics = [-2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
ztics = [750.0, 800.0, 850.0, 900.0, 950.0]
lab = ["star", "DM"]


def draw_figure(fileid, kind, Ndisk, disk_xi, disk_eta, disk_D, Nhsc, hsc_xi, hsc_eta, hsc_rad, Nnws, nws_xi, nws_eta, nws_D, nws_Derr, nws_s1_xi, nws_s1_eta, nws_s2_xi, nws_s2_eta, nws_s3_xi, nws_s3_eta, nws_s4_xi, nws_s4_eta, Ngcs, nws_gc_xi, nws_gc_eta, nws_gc_vel, nws_gc_verr):
    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".m31obs" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    time = h5file["/"].attrs["time"][0]
    useDegree = h5file["/"].attrs["useDegree"][0]

    # kind = h5file["/"].attrs["kinds"][0]
    nx = h5file["/"].attrs["nx"][0]
    ny = h5file["/"].attrs["ny"][0]
    nz = h5file["/"].attrs["nz"][0]
    nv = h5file["/"].attrs["nv"][0]

    # memory allocation for surface density maps
    xy_map = np.zeros((kind, nx, ny))
    yz_map = np.zeros((kind, ny, nz))
    xx     = np.zeros((kind, nx + 1))
    yy     = np.zeros((kind, ny + 1))
    zz     = np.zeros((kind, nz + 1))

    # memory allocation for phase-space densty map
    yv_map = np.zeros((kind, ny, nv))
    vv     = np.zeros((kind, nv + 1))

    mbh_obs = [0] * 3

    nxpanel = 3
    nypanel = 0
    Ndata = 0
    for ii in range(kind):
        if (Nskip == 0) or (ii not in skip):
            # read surface density maps
            folder = "field" + str(ii) + "/"
            xy_map[nypanel] = h5file[folder + "Sigma_xy"].value
            yz_map[nypanel] = h5file[folder + "Sigma_yz"].value
            yv_map[nypanel] = h5file[folder +     "f_yv"].value
            xx[nypanel] = h5file[folder +   "xi"].value
            yy[nypanel] = h5file[folder +  "eta"].value
            zz[nypanel] = h5file[folder +    "D"].value
            vv[nypanel] = h5file[folder + "vlos"].value

            nypanel += 1
            Ndata += 1
        if (Nmbh > 0) and (ii in wmbh):
            # read particle location
            folder = "obs" + str(ii) + "/"
            mbh_xi   = h5file[folder +   "xi"].value
            mbh_eta  = h5file[folder +  "eta"].value
            mbh_dist = h5file[folder + "dist"].value
            mbh_obs[0] = mbh_xi[0]
            mbh_obs[1] = mbh_eta[0]
            mbh_obs[2] = mbh_dist[0]
    # close the HDF5 file
    h5file.close()


    xy_map = np.maximum(xy_map, Sigma_min)
    xy_map = np.minimum(xy_map, Sigma_max)
    yz_map = np.maximum(yz_map, depth_min)
    yz_map = np.minimum(yz_map, depth_max)
    yv_map = np.maximum(yv_map, phase_min)
    yv_map = np.minimum(yv_map, phase_max)


    # plot the data
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    xmin, xmax = xx[0][0], xx[0][nx]
    ymin, ymax = yy[0][0], yy[0][ny]
    zmin, zmax = zz[0][0], zz[0][nz]
    vmin, vmax = vv[0][0], vv[0][nv]

    # adjust marker size and etcetra
    fs = fs_base / np.sqrt(nypanel)
    tl = tl_base / np.sqrt(nypanel)
    tw = tw_base / np.sqrt(nypanel)
    lw_ref = lw_ref_base / nypanel
    lw_M31disk = lw_M31disk_base / nypanel
    ms_wmbh = ms_wmbh_base / nypanel
    lw_hsc = lw_hsc_base / nypanel
    lw_nws = lw_nws_base / nypanel
    lw_gcs = lw_gcs_base / nypanel
    ms_hsc = ms_hsc_base / nypanel
    ms_nws = ms_nws_base / nypanel
    ms_gcs = ms_gcs_base / nypanel


    for jj in range(nypanel):
        img_Sigma = ax[              jj].imshow(xy_map[jj].T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = Sigma_min, vmax = Sigma_max), cmap = cmap_Sigma, aspect = "auto")
        img_depth = ax[    nypanel + jj].imshow(yz_map[jj]  , extent = [zmin, zmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = depth_min, vmax = depth_max), cmap = cmap_depth, aspect = "auto")
        img_phase = ax[2 * nypanel + jj].imshow(yv_map[jj]  , extent = [vmin, vmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = phase_min, vmax = phase_max), cmap = cmap_phase, aspect = "auto")


    for jj in range(nypanel):
        # reference circles for projected distances
        c050 = plt.Circle((0, 0), m31.kpc2degree( 50.0), ec = col_ref, lw = lw_ref, fill = False)
        c100 = plt.Circle((0, 0), m31.kpc2degree(100.0), ec = col_ref, lw = lw_ref, fill = False)
        c150 = plt.Circle((0, 0), m31.kpc2degree(150.0), ec = col_ref, lw = lw_ref, fill = False)
        c200 = plt.Circle((0, 0), m31.kpc2degree(200.0), ec = col_ref, lw = lw_ref, fill = False)
        ax[jj].add_patch(c050)
        ax[jj].add_patch(c100)
        ax[jj].add_patch(c150)
        ax[jj].add_patch(c200)

        # reference ellipse of the M31 disk
        ax[              jj].plot(disk_xi[                    ::int((Ndisk - 1) / 2)], disk_eta[                    ::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# minor axis
        ax[              jj].plot(disk_xi[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_eta[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# major axis
        ax[              jj].plot(disk_xi                                            , disk_eta                                            , "-", color = col_M31disk, linewidth = lw_M31disk)
        ax[    nypanel + jj].plot(disk_D [                    ::int((Ndisk - 1) / 2)], disk_eta[                    ::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# minor axis
        ax[    nypanel + jj].plot(disk_D [int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_eta[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], "-", color = col_M31disk, linewidth = lw_M31disk)# major axis
        ax[    nypanel + jj].plot(disk_D                                             , disk_eta                                            , "-", color = col_M31disk, linewidth = lw_M31disk)

        # reference points of observation fields
        for kk in range(Nhsc):
            circle = plt.Circle((hsc_xi[kk], hsc_eta[kk]), hsc_rad[kk], ec = col_hsc, lw = lw_hsc, fill = False)
            ax[jj].add_patch(circle)

        # distance measurements to NW stream
        ax[              jj].plot(nws_s1_xi, nws_s1_eta, "-", color = col_nws, linewidth = lw_nws)
        ax[              jj].plot(nws_s2_xi, nws_s2_eta, "-", color = col_nws, linewidth = lw_nws)
        ax[              jj].plot(nws_s3_xi, nws_s3_eta, "-", color = col_nws, linewidth = lw_nws)
        ax[              jj].plot(nws_s4_xi, nws_s4_eta, "-", color = col_nws, linewidth = lw_nws)
        ax[    nypanel + jj].plot(nws_D , nws_eta, "s", color = col_nws, markerfacecolor = "none", markersize = ms_nws)
        ax[    nypanel + jj].errorbar(nws_D, nws_eta, xerr = nws_Derr, ls = "none", ecolor = col_nws, elinewidth = lw_nws)

        # line-of-sight velocity measurements of globular clusters
        ax[              jj].plot(nws_gc_xi , nws_gc_eta, "o", color = col_gcs, markerfacecolor = "none", markersize = ms_gcs)
        ax[2 * nypanel + jj].plot(nws_gc_vel, nws_gc_eta, "o", color = col_gcs, markerfacecolor = "none", markersize = ms_gcs)
        ax[2 * nypanel + jj].errorbar(nws_gc_vel, nws_gc_eta, xerr = nws_gc_verr, ls = "none", ecolor = col_gcs, elinewidth = lw_gcs)


        # indicate wandering MBH location
        if Nmbh == 1:
            ax[              jj].plot(mbh_obs[0], mbh_obs[1], "+", color = col_wmbh, markersize = ms_wmbh)
            ax[    nypanel + jj].plot(mbh_obs[2], mbh_obs[1], "+", color = col_wmbh, markersize = ms_wmbh)



    for ii in range(nxpanel * nypanel):
        ax[ii].set_ylim([ymin, ymax])
        ax[ii].set_yticks(ytics)
        ax[ii].tick_params(axis = "both", direction = "in", color = col_grid, bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
        ax[ii].spines["bottom"].set_color(col_grid)
        ax[ii].spines[   "top"].set_color(col_grid)
        ax[ii].spines[  "left"].set_color(col_grid)
        ax[ii].spines[ "right"].set_color(col_grid)

    for jj in range(nypanel):
        ax[              jj].set_xlim([xmin, xmax])
        ax[              jj].set_xticks(xtics)
        ax[              jj].spines["left"].set_color(col_frame)
        ax[              jj].set_ylabel(r"$\eta$ ({:<})".format("degree"), fontsize = fs)
        ax[    nypanel + jj].set_xlim([zmin, zmax])
        ax[    nypanel + jj].set_xticks(ztics)
        ax[2 * nypanel + jj].set_xlim([vmin, vmax])
        # ax[2 * nypanel + jj].set_xticks(vtics)
        ax[(nxpanel - 1) * nypanel + jj].spines["right"].set_color(col_frame)

    for ii in range(nxpanel):
        ax[ii * nypanel                ].spines["bottom"].set_color(col_frame)
        ax[ii * nypanel + (nypanel - 1)].spines[   "top"].set_color(col_frame)

    ax[0          ].set_xlabel(r"$\xi$ ({:<})".format("degree"), fontsize = fs)
    ax[    nypanel].set_xlabel(r"$D$ ({:<})".format("kpc"), fontsize = fs)
    ax[2 * nypanel].set_xlabel(r"$v_\mathrm{los}$ (km s$^{-1}$)", fontsize = fs)

    for ii in range(nxpanel):
        for jj in range(nypanel):
            caption  = "(" + "{:^c}".format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ")"
            idx = ii * nypanel + jj
            if ii == 0:
                ax[idx].text(xmin + 0.03 * (xmax - xmin), ymax - 0.06 * (nypanel ** 0.5) * (ymax - ymin), caption + " " + lab[nypanel - 1 - jj], color = col_caption, fontsize = fs)
            if ii == 1:
                ax[idx].text(zmin + 0.03 * (zmax - zmin), ymax - 0.06 * (nypanel ** 0.5) * (ymax - ymin), caption, color = col_caption, fontsize = fs)
            if ii == 2:
                ax[idx].text(vmin + 0.03 * (vmax - vmin), ymax - 0.06 * (nypanel ** 0.5) * (ymax - ymin), caption, color = col_caption, fontsize = fs)



    for ii in range(nxpanel):
        # add colorbar
        x0, x1 = 1, 0
        y0, y1 = 1, 0
        # for at in ax[ii * nypanel:ii * nypanel + (nypanel - 1)]:
        for at in ax[ii * nypanel:(ii + 1) * nypanel]:
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

        colorbar_ax = fig.add_axes([x0, y1, x1 - x0, 0.05 / nypanel])

        if ii == 0:
            if useDegree:
                cbar = fig.colorbar(img_Sigma, cax = colorbar_ax, orientation = "horizontal", label = r"$\Sigma$ ($M_\odot$ deg$^{-2}$)")
            else:
                cbar = fig.colorbar(img_Sigma, cax = colorbar_ax, orientation = "horizontal", label = r"$\Sigma$ ($M_\odot$ kpc$^{-2}$)")
        if ii == 1:
            if useDegree:
                cbar = fig.colorbar(img_depth, cax = colorbar_ax, orientation = "horizontal", label = r"$\Sigma$ ($M_\odot$ deg$^{-1}$ kpc$^{-1}$)")
            else:
                cbar = fig.colorbar(img_depth, cax = colorbar_ax, orientation = "horizontal", label = r"$\Sigma$ ($M_\odot$ kpc$^{-2}$)")
        if ii == 2:
            if useDegree:
                cbar = fig.colorbar(img_phase, cax = colorbar_ax, orientation = "horizontal", label = r"$f$ ($M_\odot$ deg$^{-1}$ km$^{-1}$ s)")
            else:
                cbar = fig.colorbar(img_phase, cax = colorbar_ax, orientation = "horizontal", label = r"$f$ ($M_\odot$ kpc$^{-1}$ km$^{-1}$ s)")
        cbar.solids.set_edgecolor("face")
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), fontsize=fs)
        cbar.ax.set_xlabel(cbar.ax.get_xlabel(), fontsize=fs)
        colorbar_ax.get_xaxis().set_ticks_position("top")
        colorbar_ax.get_xaxis().set_label_position("top")


    # add current time
    if not monochrome:
        if nypanel == 1:
            fig.suptitle(r"$t = {:.3f}$ {:<}".format(time / 1000, "Gyr"), x = 0.37, y = 1.0, fontsize = fs)
        else:
            fig.suptitle(r"$t = {:.3f}$ {:<}".format(time / 1000, "Gyr"), y = 1.0, fontsize = fs)


    # save figures
    figname = "fig/" + filename + "_map"
    if monochrome:
        figname += "_mono"
    figname += snapshot
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

# set reference points of observations by Komiyama et al. (2018)
Nhsc, hsc_xi, hsc_eta, hsc_rad, Nnws, nws_xi, nws_eta, nws_D, nws_Dep, nws_Dem, nws_s1_xi, nws_s1_eta, nws_s2_xi, nws_s2_eta, nws_s3_xi, nws_s3_eta, nws_s4_xi, nws_s4_eta = m31.NWS_Komiyama2018()
nws_Derr = np.zeros((2, Nnws))
for ii in range(Nnws):
    nws_Derr[0][ii] = nws_Dem[ii]
    nws_Derr[1][ii] = nws_Dep[ii]

# set reference points of line-of-sight velocity measurements of GCs by Veljanoski et al. (2014)
Ngcs, nws_gc_xi, nws_gc_eta, nws_gc_vel, nws_gc_vep, nws_gc_vem = m31.NWS_velocity()
nws_gc_verr = np.zeros((2, Ngcs))
for ii in range(Ngcs):
    nws_gc_verr[0][ii] = nws_gc_vem[ii]
    nws_gc_verr[1][ii] = nws_gc_vep[ii]

cmap_Sigma = utils.generate_cmap(["darkblue", "deepskyblue", "lime", "yellow", "red", "magenta", "white"])
cmap_depth = "jet"
cmap_phase = "hot"
if monochrome:
    cmap_Sigma = "gray_r"
    cmap_depth = "gray_r"
    cmap_phase = "gray_r"
cores = int(np.ceil(mp.cpu_count() / 2))
pool = mp.Pool(cores)
args = [(ii, Nkind, Ndisk, disk_xi, disk_eta, disk_D, Nhsc, hsc_xi, hsc_eta, hsc_rad, Nnws, nws_xi, nws_eta, nws_D, nws_Derr, nws_s1_xi, nws_s1_eta, nws_s2_xi, nws_s2_eta, nws_s3_xi, nws_s3_eta, nws_s4_xi, nws_s4_eta, Ngcs, nws_gc_xi, nws_gc_eta, nws_gc_vel, nws_gc_verr) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
