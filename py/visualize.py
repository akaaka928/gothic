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


outputPDF = False


# filename = "m31"
# Nskip = 1
# skip = [0]
# init = 0
# last = 47
# fmin, fmax = 1.0e+7, 1.0e+10
# xmin, xmax = -15.0, 15.0
# ymin, ymax = -15.0, 15.0
# zmin, zmax = -15.0, 15.0
# xtics = [-10, -5, 0, 5, 10]
# ytics = [-10, -5, 0, 5, 10]

filename = "k17disk"
Nskip = 2
skip = [0, 2]
init = 0
last = 399
fmin, fmax = 1.0e-7, 1.0e-1
lab = ["halo", "bulge", "MBH", "disk"]

# filename = "halocusp_run"
# Nskip = 1
# skip = [0]
# init = 0
# last = 140
# # last = 0
# fmin, fmax = 1.0e-6, 1.0e-1
# lab = ["halo", "core", "GC"]

# filename = "hernquist"
# Nskip = 0
# init = 0
# last = 47
# # last = 0
# fmin, fmax = 1.0e-3, 1.0e+1
# lab = ["halo"]



pt = ["o", "s", "^", "D"]
ls = ["-", ":", "-.", "--"]
col = ["black", "red", "blue", "magenta"]



def draw_figure(fileid, kind):
    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".plt" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    length_unit = h5file["/"].attrs["length_astro_unit_name"][0]
    time_unit = h5file["/"].attrs["time_astro_unit_name"][0]
    mass_unit = h5file["/"].attrs["mass_astro_unit_name"][0]
    velocity_unit = h5file["/"].attrs["velocity_astro_unit_name"][0]
    density_unit = h5file["/"].attrs["density_astro_unit_name"][0]
    col_density_unit = h5file["/"].attrs["col_density_astro_unit_name"][0]
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


    fig_rho = utils.set_figure(1, 1)
    fig_enc = utils.set_figure(1, 1)
    fig_sig = utils.set_figure(1, 1)
    fig_Sig = utils.set_figure(1, 1)
    fig_ver = utils.set_figure(1, 1)
    fig_sig2D = utils.set_figure(1, 1)
    ax_rho = [0]
    ax_enc = [0]
    ax_sig = [0]
    ax_Sig = [0]
    ax_ver = [0]
    ax_sig2D = [0]
    utils.locate_panels(fig_rho, ax_rho, 1, 1, True, True)
    utils.locate_panels(fig_enc, ax_enc, 1, 1, True, True)
    utils.locate_panels(fig_sig, ax_sig, 1, 1, True, True)
    utils.locate_panels(fig_Sig, ax_Sig, 1, 1, True, True)
    utils.locate_panels(fig_ver, ax_ver, 1, 1, True, True)
    utils.locate_panels(fig_sig2D, ax_sig2D, 1, 1, True, True)


    nxpanel = 0
    nypanel = 2
    Ndata = 0
    for ii in range(kind):
        if (Nskip == 0) or (ii not in skip):
            num = h5file["attr" + str(ii)].attrs["number"][0]
            if (num > 1):
                # read surface density maps
                folder = "field" + str(ii) + "/"
                xy_map[nxpanel] = h5file[folder + "Sigma_xy"].value
                xz_map[nxpanel] = h5file[folder + "Sigma_zx"].value.T
                xx[nxpanel] = h5file[folder + "x"].value
                yy[nxpanel] = h5file[folder + "y"].value
                zz[nxpanel] = h5file[folder + "z"].value

                nxpanel += 1
                Ndata += 1

            folder = "rad" + str(ii) + "/"
            rad = h5file[folder + "rad"].value
            rho = h5file[folder + "rho"].value
            enc = h5file[folder + "enc"].value
            sig = h5file[folder + "sig"].value

            folder = "hor" + str(ii) + "/"
            hor = h5file[folder + "hor"].value
            ver = h5file[folder + "height"].value
            Sig = h5file[folder + "Sigma"].value
            sigR = h5file[folder + "sigR"].value
            sigp = h5file[folder + "sigp"].value
            sigz = h5file[folder + "sigz"].value

            # measure of the thickness, according to Rodionov & Sotnikova (2013)
            ver *= 1.82

            ax_rho[0].plot(rad, rho, pt[ii], linestyle = ls[ii], color = col[ii], label = lab[ii])
            ax_enc[0].plot(rad, enc, pt[ii], linestyle = ls[ii], color = col[ii], label = lab[ii])
            ax_sig[0].plot(rad, sig, pt[ii], linestyle = ls[ii], color = col[ii], label = lab[ii])
            ax_Sig[0].plot(hor, Sig, pt[ii], linestyle = ls[ii], color = col[ii], label = lab[ii])
            ax_ver[0].plot(hor, ver, pt[ii], linestyle = ls[ii], color = col[ii], label = lab[ii])
            ax_sig2D[0].plot(hor, sigR, pt[0], linestyle = ls[ii], color = col[ii], label = r"$\sigma_R$" + " (" + lab[ii] + ")")
            ax_sig2D[0].plot(hor, sigp, pt[1], linestyle = ls[ii], color = col[ii], label = r"$\sigma_p$" + " (" + lab[ii] + ")")
            ax_sig2D[0].plot(hor, sigz, pt[2], linestyle = ls[ii], color = col[ii], label = r"$\sigma_z$" + " (" + lab[ii] + ")")




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
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    xmin, xmax = xx[0][0], xx[0][nx]
    ymin, ymax = yy[0][0], yy[0][ny]
    zmin, zmax = zz[0][0], zz[0][nz]

    head = 0
    if summary:
        img = ax[1].imshow(xz_map[nxpanel - 1].T, extent = [xmin, xmax, zmin, zmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img = ax[0].imshow(xy_map[nxpanel - 1].T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        head = 1

    for ii in range(Ndata):
        img = ax[(head + ii) * nypanel + 1].imshow(xz_map[ii].T, extent = [xmin, xmax, zmin, zmax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        img = ax[(head + ii) * nypanel    ].imshow(xy_map[ii].T, extent = [xmin, xmax, ymin, ymax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")


    for ii in range(nxpanel * nypanel):
        ax[ii].set_xlim([xmin, xmax])
        # ax[ii].set_xticks(xtics)
        # ax[ii].set_yticks(ytics)
        ax[ii].tick_params(axis = "both", direction = "in", color = "white", bottom = True, top = True, left = True, right = True)
        ax[ii].spines["bottom"].set_color("white")
        ax[ii].spines[   "top"].set_color("white")
        ax[ii].spines[  "left"].set_color("white")
        ax[ii].spines[ "right"].set_color("white")

    for ii in range(nxpanel):
        ax[ii * nypanel + 1].set_ylim([zmin, zmax])
        ax[ii * nypanel    ].set_ylim([ymin, ymax])
        ax[ii * nypanel    ].set_xlabel(r"$x$ ({:<})".format(length_unit.decode("UTF-8")))
        ax[ii * nypanel + 1].spines[   "top"].set_color("black")
        ax[ii * nypanel    ].spines["bottom"].set_color("black")


    ax[0].set_ylabel(r"$y$ ({:<})".format(length_unit.decode("UTF-8")))
    ax[1].set_ylabel(r"$z$ ({:<})".format(length_unit.decode("UTF-8")))
    ax[0].spines["left"].set_color("black")
    ax[1].spines["left"].set_color("black")
    ax[(nxpanel - 1) * nypanel + 1].spines[ "right"].set_color("black")
    ax[(nxpanel - 1) * nypanel    ].spines[ "right"].set_color("black")



    ax_rho[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_enc[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_sig[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_Sig[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_ver[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_sig2D[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode("UTF-8")))

    ax_rho[0].set_ylabel(r"$\rho$ ({:<})".format(density_unit.decode("UTF-8")))
    ax_enc[0].set_ylabel(r"$M_\mathrm{enc}$" + r"({:<})".format(mass_unit.decode("UTF-8")))
    ax_sig[0].set_ylabel(r"$\sigma_r$ ({:<})".format(velocity_unit.decode("UTF-8")))
    ax_Sig[0].set_ylabel(r"$\Sigma$ ({:<})".format(col_density_unit.decode("UTF-8")))
    ax_ver[0].set_ylabel(r"$1.82\, \mathrm{median}(|z|)$" + r"({:<})".format(length_unit.decode("UTF-8")))
    ax_sig2D[0].set_ylabel(r"$\sigma$ ({:<})".format(velocity_unit.decode("UTF-8")))

    ax_rho[0].loglog()
    ax_enc[0].loglog()
    ax_Sig[0].loglog()
    ax_sig[0].semilogx()
    ax_ver[0].semilogx()
    ax_sig2D[0].semilogx()

    ax_rho[0].grid()
    ax_enc[0].grid()
    ax_sig[0].grid()
    ax_Sig[0].grid()
    ax_ver[0].grid()
    ax_sig2D[0].grid()


    # add legends
    handles, labels = ax_rho[0].get_legend_handles_labels()
    # ax_rho[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best')
    ax_rho[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')
    handles, labels = ax_enc[0].get_legend_handles_labels()
    ax_enc[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')
    handles, labels = ax_sig[0].get_legend_handles_labels()
    ax_sig[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')
    handles, labels = ax_Sig[0].get_legend_handles_labels()
    ax_Sig[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')
    handles, labels = ax_ver[0].get_legend_handles_labels()
    ax_ver[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')
    handles, labels = ax_sig2D[0].get_legend_handles_labels()
    ax_sig2D[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')


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
    cbar = fig.colorbar(img, cax = colorbar_ax, label = r"$\Sigma$ ({:<})".format(col_density_unit.decode("UTF-8")))
    # cbar = fig.colorbar(img, cax = colorbar_ax, label = r"$\Sigma$ ($M_\odot$ kpc$^{-2}$)")
    cbar.solids.set_edgecolor("face")

    # add current time
    fig.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    # fig.suptitle(r"$t = {:.2f}$ {:<}".format(time / 1000, "Gyr"))

    fig_rho.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    fig_enc.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    fig_sig.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    fig_Sig.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    fig_ver.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    fig_sig2D.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))



    # save figures
    figname = "fig/" + filename + "_map" + snapshot
    fig.savefig(figname + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_rho.savefig("fig/" + filename + "_rho" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_enc.savefig("fig/" + filename + "_enc" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_sig.savefig("fig/" + filename + "_sig" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_Sig.savefig("fig/" + filename + "_Sig" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_ver.savefig("fig/" + filename + "_ver" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_sig2D.savefig("fig/" + filename + "_sigRz" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if outputPDF:
        fig.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_rho.savefig("fig/" + filename + "_rho" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_enc.savefig("fig/" + filename + "_enc" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_sig.savefig("fig/" + filename + "_sig" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_Sig.savefig("fig/" + filename + "_Sig" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_ver.savefig("fig/" + filename + "_ver" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_sig2D.savefig("fig/" + filename + "_sigRz" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.close("all")



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
txtfile = open("doc/" + filename + ".summary.txt", "r")
unit = txtfile.readline()
tmp = txtfile.readline()
Nkind = int(tmp[0])
# Nsphe = int(tmp[2])
txtfile.close()

my_cmap = utils.generate_cmap(["darkblue", "deepskyblue", "lime", "yellow", "red", "magenta", "white"])
# cores = mp.cpu_count()
cores = int(np.ceil(mp.cpu_count() / 2))
pool = mp.Pool(cores)
args = [(ii, Nkind) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
