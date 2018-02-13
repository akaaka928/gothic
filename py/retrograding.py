# show retrograding fraction f as a function of R
import numpy as np
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import os.path as path
import multiprocessing as mp

import utils as utils


# specify plot target
# filename = "m31"
filename = "sat"
# init = 0
# last = 47
Ncrit = 2048 # number of particles for estimating physical quantities
Nmin = 16 # minimum data points for visualization

tag = ["disk"]
pt = ["o", "s", "^", "D"]
ls = ["-", ":", "-.", "--"]
col = ["black", "red", "blue", "magenta"]

# set plot range
Rmin, Rmax = 0.0, 25.0


# set number of panels
nxpanel, nypanel = 1, 1


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 28


# read fundamental property of the system
target = "doc/" + filename + ".summary.txt"
txtfile = open(target, "r")
unit = txtfile.readline()
# print(unit)
# Nkind, Nsphe = txtfile.readline()
tmp = txtfile.readline()
Nkind = int(tmp[0])
Nsphe = int(tmp[2])
# print(tmp)
# print(Nkind)
# print(Nsphe)
num = [0] * Nkind
Ntot = 0
for ii in range(Nkind):
    tmp = txtfile.readline()
    num[ii] = int(tmp)
    Ntot += num[ii]
txtfile.close()
# print(num)

head = [0] * Nkind
head[0] = 0
for ii in range(Nkind - 1):
    head[ii + 1] = head[ii] + num[ii]

Ndisk = Nkind - Nsphe
Ndat_disk = [0] * Ndisk
head_disk = [0] * Ndisk
head_disk[0] = 0
Nmax = int(np.ceil((Ntot / Ncrit) + Ndisk * Nmin))
RR_plot = [0] * Nmax
fr_plot = [0] * Nmax


fig = utils.set_figure(nxpanel, nypanel)
ax = [0] * nxpanel * nypanel
utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)


# read fundamental property of disk component(s)
target = "dat/" + filename + ".tmp0.h5"
if path.isfile(target) == True:
    data_file = h5py.File(target, "r")

    # read particle position and velocity
    folder = "nbody/"
    pos = data_file[folder + "position"].value
    vel = data_file[folder + "velocity"].value
    idx = data_file[folder + "index"].value

    Ntot = data_file["nbody/"].attrs["number"][0]
    posx = np.empty(Ntot)
    posy = np.empty(Ntot)
    velx = np.empty(Ntot)
    vely = np.empty(Ntot)
    for ii in range(Ntot):
        posx[ii] = pos[4 * ii    ]
        posy[ii] = pos[4 * ii + 1]
        velx[ii] = vel[4 * ii    ]
        vely[ii] = vel[4 * ii + 1]

    R2 = posx * posx + posy * posy
    Ri = 1.0 / np.sqrt(1.0e-100 + R2)
    RR = R2 * Ri
    lz = posx * vely - posy * velx

    src = np.argsort(idx)

    for ii in range(Ndisk):
        kk = Nsphe + ii

        RR_disk = np.empty(num[kk])
        lz_disk = np.empty(num[kk])
        for jj in range(num[kk]):
            ll = src[head[kk] + jj]
            RR_disk[jj] = RR[ll]
            lz_disk[jj] = lz[ll]

        src_disk = np.argsort(RR_disk)

        nunit = Ncrit
        if num[kk] < (Ncrit * Nmin):
            nunit = int(np.floor(num[kk] / Nmin))
        Ndat_disk[ii] = int(np.floor(num[kk] / nunit))

        RR_sort = np.empty(nunit)
        lz_sort = np.empty(nunit)

        for jj in range(Ndat_disk[ii]):
            Nneg = 0
            Ntot = 0
            for ll in range(nunit):
                mm = src_disk[jj * nunit + ll]
                RR_sort[ll] = RR_disk[mm]
                lz_sort[ll] = lz_disk[mm]
                Ntot += 1
                if lz_sort[ll] < 0:
                    Nneg += 1

            RR_plot[head_disk[ii] + jj] = np.median(RR_sort)
            fr_plot[head_disk[ii] + jj] = Nneg / Ntot

        if ii != (Ndisk - 1):
            head_disk[ii + 1] = head_disk[ii] + Ndat_disk[ii]

    data_file.close()


    # plot Lz profile of disk component(s)
    for kk in range(Ndisk - 1, -1, -1):
        ax[0].plot(RR_plot[head_disk[kk] : head_disk[kk] + Ndat_disk[kk]], fr_plot[head_disk[kk] : head_disk[kk] + Ndat_disk[kk]], ls[kk], color = col[kk], label = tag[kk])
        # ax[0].plot(RR_plot[head_disk[kk] : head_disk[kk] + Ndat_disk[kk]], fr_plot[head_disk[kk] : head_disk[kk] + Ndat_disk[kk]], pt[kk], color = col[kk], label = tag[kk], markerfacecolor = "none")

    # set plot range
    # ax[0].set_xlim([Rmin, Rmax])
    # ax[0].set_ylim([Lmin, Lmax])
    # ax[0].semilogy()
    ax[0].semilogx()
    # ax[0].loglog()
    ax[0].grid()

    # set label
    ax[0].set_xlabel(r"$R$ (kpc)")
    ax[0].set_ylabel(r"Retrograding fraction")

    # # set legend
    # handles, labels = ax[0].get_legend_handles_labels()
    # # ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, prop = {'size' : 14}, loc = 'best')
    # ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best')

    # output figure
    figname = "fig/" + filename + "_retrograding"
    plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")


    # clear figure
    ax[0].cla()






else:
    print("{:<} was not found.".format(target))
    sys.exit()



