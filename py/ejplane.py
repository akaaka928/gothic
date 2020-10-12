import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm # for logarithmic plot in imshow
import matplotlib.ticker as ticker

import os.path as path
import multiprocessing as mp

import utils as utils


outputPDF = False

closeUp = False


fs_base = 24
tl_base = 6.0
tw_base = 1.0


senergy2astro = 1.0e+10


# filename = "m15ra1_4_run"
# init = 0
# last = 80
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["DM halo", "bulge"]
# fmin, fmax = 3.1e+0 / senergy2astro, 3.1e+3 / senergy2astro

# filename = "w3iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w3ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w3ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.5e+5 / senergy2astro

# filename = "w3ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 2.0e+5 / senergy2astro

# filename = "w3ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "w3ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+4 / senergy2astro, 1.0e+6 / senergy2astro

# filename = "w5iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w5ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w5ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w5ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 2.0e+5 / senergy2astro

# filename = "w5ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "w5ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 8.0e+5 / senergy2astro

# filename = "w7iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w7ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w7ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w7ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "w7ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "w7ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+6 / senergy2astro

# filename = "m12iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 4.0e+3 / senergy2astro

# filename = "m12ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 5.0e+3 / senergy2astro

filename = "m12ra1"
init = 0
last = 140
# last = 0
Nskip = 0
skip = [0]
lab = ["halo"]
fmin, fmax = 1.0e+2 / senergy2astro, 1.0e+4 / senergy2astro

# filename = "m12ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 1.0e+4 / senergy2astro

# filename = "m12ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 3.1e+4 / senergy2astro

# filename = "m12ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+0 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "m09iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "m09ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "m09ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "m09ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "m09ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "m09ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "a00iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 3.1e+4 / senergy2astro

# filename = "a00ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 3.1e+4 / senergy2astro

# filename = "a00ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 3.1e+4 / senergy2astro

# filename = "a00ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "a00ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 3.1e+5 / senergy2astro


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", "--", ":", "-."]
Ncol, col = utils.set_color_palette_for_color_universal_design()


def draw_figure(fileid, kind, sphe):
    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".plt" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    time = h5file["/"].attrs["time"][0]

    nE = h5file["/"].attrs["nE"][0]
    nJ = h5file["/"].attrs["nJ"][0]

    # memory allocation for EJ-maps
    EJful = np.zeros((kind + 1, nJ, nE))
    EJsml = np.zeros((kind + 1, nJ, nE))
    EE    = np.zeros((kind + 1, nE + 1))
    Jful  = np.zeros((kind + 1, nJ + 1))
    Jsml  = np.zeros((kind + 1, nJ + 1))


    nxpanel = 1
    if closeUp:
        nxpanel = 2
    nypanel = 0
    Ndata = 0
    for ii in range(kind):
        num = h5file["attr" + str(ii)].attrs["number"][0]
        if (num > 1):
            if (Nskip == 0) or (ii not in skip):
                # EJ-maps
                folder = "EJplane" + str(ii) + "/"
                EJful[nypanel] = h5file[folder + "EJful"] / senergy2astro
                EJsml[nypanel] = h5file[folder + "EJsml"] / senergy2astro
                EE[nypanel] = h5file[folder + "E"] * senergy2astro
                Jful[nypanel] = h5file[folder + "Jful"]
                Jsml[nypanel] = h5file[folder + "Jsml"]

                nypanel += 1
                Ndata += 1

    # close the HDF5 file
    h5file.close()

    summary = False
    if nypanel > 1:
        for ii in range(nypanel):
            EJful[nypanel] += EJful[ii]
            EJsml[nypanel] += EJsml[ii]

        nypanel += 1
        summary = True


    EJful = np.maximum(EJful, fmin)
    EJful = np.minimum(EJful, fmax)
    EJsml = np.maximum(EJsml, fmin)
    EJsml = np.minimum(EJsml, fmax)


    # plot the data
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    Emin, Emax = EE[0][0], EE[0][nE]
    Jfulmin, Jfulmax = Jful[0][0], Jful[0][nJ]
    Jsmlmin, Jsmlmax = Jsml[0][0], Jsml[0][nJ]

    # adjust marker size and etcetra
    fs = fs_base / np.sqrt(max([nxpanel, nypanel]))
    tl = tl_base / np.sqrt(max([nxpanel, nypanel]))
    tw = tw_base / np.sqrt(max([nxpanel, nypanel]))

    if summary:
        img = ax[    nypanel - 1].imshow(EJful[nypanel - 1].T, extent = [Jfulmin, Jfulmax, Emin, Emax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        if closeUp:
            img = ax[2 * nypanel - 1].imshow(EJsml[nypanel - 1].T, extent = [Jsmlmin, Jsmlmax, Emin, Emax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

    for ii in range(Ndata):
        img = ax[ii          ].imshow(EJful[ii].T, extent = [Jfulmin, Jfulmax, Emin, Emax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        if closeUp:
            img = ax[ii + nypanel].imshow(EJsml[ii].T, extent = [Jsmlmin, Jsmlmax, Emin, Emax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

    for ii in range(nxpanel * nypanel):
        ax[ii].set_ylim([Emin, Emax])
        ax[ii].tick_params(axis = "both", direction = "in", color = "white", bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
        ax[ii].spines["bottom"].set_color("white")
        ax[ii].spines[   "top"].set_color("white")
        ax[ii].spines[  "left"].set_color("white")
        ax[ii].spines[ "right"].set_color("white")
        ax[ii].yaxis.set_major_formatter(ticker.FuncFormatter(utils.scientific))

    for ii in range(nypanel):
        ax[ii          ].set_xlim([Jfulmin, Jfulmax])
        ax[ii          ].spines[ "left"].set_color("black")
        ax[ii          ].set_ylabel(r"$E$~(\si{erg.g^{-1}})", fontsize = fs)
        if closeUp:
            ax[ii + nypanel].set_xlim([Jsmlmin, Jsmlmax])
            ax[ii + nypanel].spines["right"].set_color("black")

    for ii in range(nxpanel):
        ax[(ii + 1) * nypanel - 1].spines["top"].set_color("black")
        ax[ii * nypanel].spines["bottom"].set_color("black")
        ax[ii * nypanel].set_xlabel(r"$J$~(\si{kpc.km.s^{-1}})", fontsize = fs)


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
    cbar = fig.colorbar(img, cax = colorbar_ax, format=ticker.FuncFormatter(utils.scientific), label = r"$f$ (" + r"\si{M_\odot.kpc^{-1}.km^{-1}.s.erg^{-1}.g}" + r")")
    cbar.solids.set_edgecolor("face")

    # add current time
    fig.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})", fontsize = fs)


    # save figures
    figname = "fig/" + filename + "_EJmap" + snapshot
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
cores = int(np.ceil(mp.cpu_count() / 2))
pool = mp.Pool(cores)
args = [(ii, Nkind, Nsphe) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
