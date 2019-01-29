import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import os.path as path
import multiprocessing as mp

import utils as utils


outputPDF = False
senergy2astro = 1.0e+10


filename = "m15ra1_4_run"
init = 0
last = 80
# last = 0
lab = ["DM halo", "bulge"]

nxpanel, nypanel = 1, 1


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", "--", ":", "-."]
Ncol, col = utils.set_color_palette_for_color_universal_design()


def draw_figure(fileid, kind, Emin, Emax, Jmin, Jmax):
    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".plt" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    time = h5file["/"].attrs["time"][0]

    nypanel = kind
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)


    for ii in range(kind):
        folder = "particle" + str(ii) + "/"
        Etot = h5file[folder + "Etot"].value * senergy2astro
        Jtot = h5file[folder + "Jtot"].value
        ax[ii * nypanel].plot(Jtot, Etot, ",", color = "black", rasterized = True)

    # close the HDF5 file
    h5file.close()

    ax[0].set_xlabel(r"$J$~(\si{kpc.km.s^{-1}})")

    for ii in range(kind):
        ax[ii * nypanel].set_ylabel(r"$E$~(\si{erg.g^{-1}})")

    if Jmin <= 0.0:
        Jmin = 1.0e-7 * Jmax
    for ii in range(nxpanel * nypanel):
        ax[ii].set_xlim(utils.scale_axis(Jmin, Jmax, True))
        ax[ii].set_ylim(utils.scale_axis(senergy2astro * Emin, senergy2astro * Emax, False))
        ax[ii].semilogx()
        ax[ii].grid()

    # set caption
    for ii in range(nxpanel):
        for jj in range(nypanel):
            idx = jj + ii * nxpanel
            caption  = "(" + "{:^c}".format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ")"
            caption += " " + lab[jj]
            ax[idx].text(0.4 * Jmin, 0.6 * senergy2astro * Emax, caption)

    # add current time
    fig.suptitle(r"$t = {:.3f}$".format(time) + r"~(\si{Myr})")

    # save figures
    fig.savefig("fig/" + filename + "_losscone" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if outputPDF:
        fig.savefig("fig/" + filename + "_losscone" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.close("all")


def wrapper(argv):
    return draw_figure(*argv)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rcParams['text.latex.preamble'] = [r"\usepackage{siunitx}"]

# set font size
plt.rcParams['font.size'] = 18
# plt.rcParams['font.size'] = 16

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


# read minimum and maximum of each quantity
txtfile = open("dat/" + filename + ".minmax.txt", "r")
line = txtfile.readline()
item = line.split("\t")
radmin, radmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
rhomin, rhomax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
encmin, encmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
sigmin, sigmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
betmin, betmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
hormin, hormax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
Sigmin, Sigmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
sigRmin, sigRmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
sigpmin, sigpmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
sigzmin, sigzmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
diskmin, diskmax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
Emin, Emax = float(item[0]), float(item[1])
line = txtfile.readline()
item = line.split("\t")
Jmin, Jmax = float(item[0]), float(item[1])
txtfile.close()


# read number of all component(s)
txtfile = open("doc/" + filename + ".summary.txt", "r")
unit = int(txtfile.readline())
line = txtfile.readline()
item = line.split("\t")
Nkind = int(item[0])
Nsphe = int(item[1])
txtfile.close()


cores = int(np.ceil(mp.cpu_count() / 2))
pool = mp.Pool(cores)
args = [(ii, Nkind - 1, Emin, Emax, Jmin, Jmax) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
