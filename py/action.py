# show angle-action variables for spherical potentials (Jr, Jtheta, Jphi), these are useful for axisymmetric disks
# plot J as a function of R
import numpy as np
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import os.path as path
import multiprocessing as mp

import utils as utils


# specify plot target
# filename = "cb17"
filename = "m31"
init = 0
last = 47
# last = 0

# tag = ["dark matter halo", "bulge", "thick disc", "thin disc"]
tag = ["dark matter halo", "bulge", "disk"]
pt = ["D", "o", "s", "^"]
ls = ["--", "-", ":", "-."]
col = ["magenta", "black", "red", "blue"]

# set plot range
# Rmin, Rmax = 0.0, 25.0
Rmin, Rmax = 1.0e-1, 25.0


# set number of panels
nxpanel, nypanel = 1, 3




def draw_figures(fileid, head, tail, xx0, Jr0, Jt0, Jp0):
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    # for fileid in range(init, last + 1, 1):
    # get snapshot ID
    snapshot = "{:03d}".format(fileid)


    # obtain angle-action variables
    h5file = h5py.File("dat/" + filename + ".action" + snapshot + ".h5", "r")
    kinds = h5file["/"].attrs["kinds"][0]
    skind = h5file["/"].attrs["skinds"][0]
    nskip = h5file["/"].attrs["skipped kinds"][0]
    time = h5file["/"].attrs["time"][0]
    time_unit = h5file["/"].attrs["time_astro_unit_name"][0]
    length_unit = h5file["/"].attrs["length_astro_unit_name"][0]


    for kk in range(kinds - 1, nskip - 1, -1):
        ax[0].plot(xx0[head[kk] : tail[kk]], Jp0[head[kk] : tail[kk]], ls[kk], color = col[kk], label = tag[kk] + " (" + r"$t = 0$~{:<}".format(time_unit.decode('UTF-8')) + ")")
        ax[1].plot(xx0[head[kk] : tail[kk]], Jt0[head[kk] : tail[kk]], ls[kk], color = col[kk], label = tag[kk] + " (" + r"$t = 0$~{:<}".format(time_unit.decode('UTF-8')) + ")")
        ax[2].plot(xx0[head[kk] : tail[kk]], Jr0[head[kk] : tail[kk]], ls[kk], color = col[kk], label = tag[kk] + " (" + r"$t = 0$~{:<}".format(time_unit.decode('UTF-8')) + ")")


    for kk in range(kinds - 1, nskip - 1, -1):
        folder = "data" + str(kk) + "/"

        Jr = h5file[folder + "Jr"    ].value
        Jt = h5file[folder + "Jtheta"].value
        Jp = h5file[folder + "Jphi"  ].value

        if kk < skind:
            xx = h5file[folder + "r"].value
        else:
            xx = h5file[folder + "R"].value

        ax[0].plot(xx, Jp, pt[kk], color = col[kk], label = tag[kk] + " (" + r"$t = {:.0f}$~{:<}".format(time, time_unit.decode('UTF-8')) + ")", markerfacecolor = "none")
        ax[1].plot(xx, Jt, pt[kk], color = col[kk], label = tag[kk] + " (" + r"$t = {:.0f}$~{:<}".format(time, time_unit.decode('UTF-8')) + ")", markerfacecolor = "none")
        ax[2].plot(xx, Jr, pt[kk], color = col[kk], label = tag[kk] + " (" + r"$t = {:.0f}$~{:<}".format(time, time_unit.decode('UTF-8')) + ")", markerfacecolor = "none")

    h5file.close()


    ax[0].set_ylim([-500.0, 3990.0])
    ax[1].set_ylim([1.0, 300.0])
    ax[2].set_ylim([1.0, 300.0])

    ax[0].semilogx()
    ax[1].loglog()
    ax[2].loglog()

    for ii in range(nxpanel * nypanel):
        # set plot range
        ax[ii].set_xlim([Rmin, Rmax])
        # ax[ii].set_ylim([Jmin, Jmax])
        # ax[ii].semilogy()
        # ax[ii].semilogx()
        # ax[ii].loglog()
        ax[ii].grid()

    # set label
    ax[0].set_xlabel(r"$r$ or $R$ ({:<})".format(length_unit.decode('UTF-8')))

    ax[0].set_ylabel(r"$J_\phi$ ({:<}".format(length_unit.decode('UTF-8')) + r" km s$^{-1}$)")
    ax[1].set_ylabel(r"$J_\theta$ ({:<}".format(length_unit.decode('UTF-8')) + r" km s$^{-1}$)")
    ax[2].set_ylabel(r"$J_r$ ({:<}".format(length_unit.decode('UTF-8')) + r" km s$^{-1}$)")

    # set legend
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = "best")

    # # fig.suptitle(r"$t = {:.0f}$ {:<}".format(time, time_unit.decode('UTF-8')))
    # ax[nypanel - 1].set_title(r"$t = {:.0f}$ {:<}".format(time, time_unit.decode('UTF-8')))

    # output figure
    figname = "fig/" + filename + "_action_" + snapshot
    plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")


    # clear figure
    for ii in range(nxpanel * nypanel):
        ax[ii].cla()


def wrapper(argv):
    return draw_figures(*argv)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 28

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


h5file = h5py.File("dat/" + filename + ".action000.h5", "r")
kinds = h5file["/"].attrs["kinds"][0]
skind = h5file["/"].attrs["skinds"][0]
nskip = h5file["/"].attrs["skipped kinds"][0]
head = [0] * kinds
tail = [0] * kinds
head[0] = 0
ndat = 0
for kk in range(nskip):
    tail[kk] = 0
    head[kk + 1] = 0
for kk in range(nskip, kinds, 1):
    num = h5file["/data" + str(kk) + "/"].attrs["number"][0]
    ndat += num
    tail[kk] = head[kk] + num
    if kk + 1 < kinds:
        head[kk + 1] = tail[kk]
xx0 = [0] * ndat
Jr0 = [0] * ndat
Jt0 = [0] * ndat
Jp0 = [0] * ndat

for kk in range(nskip, kinds, 1):
    folder = "data" + str(kk) + "/"
    num = h5file[folder].attrs["number"][0]

    Jr = h5file[folder + "Jr"    ].value
    Jt = h5file[folder + "Jtheta"].value
    Jp = h5file[folder + "Jphi"  ].value

    if kk < skind:
        xx = h5file[folder + "r"].value
    else:
        xx = h5file[folder + "R"].value

    for ii in range(num):
        xx0[head[kk] + ii] = xx[ii]
        Jr0[head[kk] + ii] = Jr[ii]
        Jt0[head[kk] + ii] = Jt[ii]
        Jp0[head[kk] + ii] = Jp[ii]
h5file.close()


cores = mp.cpu_count()
# cores = int(np.ceil(cores / 2))

pool = mp.Pool(cores)
args = [(ii, head, tail, xx0, Jr0, Jt0, Jp0) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
