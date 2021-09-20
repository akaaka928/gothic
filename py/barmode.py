# show t = T / |W|, where T is rotational kinetic energy and W is potential energy
# t > ~0.14 is unstable to a bar-like mode instability
# plot t as a function of time
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
Ncrit = 2048 # number of particles for estimating physical quantities
Nmin = 16 # minimum data points for visualization

# tag = ["total", "thick disc", "thin disc"]
tag = ["total", "disk"]
pt = ["o", "s", "^", "D"]
ls = ["-", ":", "-.", "--"]
col = ["black", "red", "blue", "magenta"]

# set plot range
# Rmin, Rmax = 0.0, 25.0
Rmin, Rmax = 1.0e-1, 25.0
Lmin, Lmax = 1.0e+5, 1.0e+9
tmin, tmax = 0.0, 0.5


# set number of panels
nxpanel, nypanel = 1, 1




def draw_figures(fileid, Nkind, Nsphe, Ndisk, mass_unit, length_unit, velocity_unit, time_unit, head0, Ndat0, RR0, Lz0):
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    # for fileid in range(init, last + 1, 1):
    # get snapshot ID
    snapshot = "{:03d}".format(fileid)

    # get file name
    input_file = "dat/" + filename + ".split" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    Ntot = h5file["/"].attrs["number"][0]
    time = h5file["/"].attrs["time"][0]


    # memory allocation for angular momentum profile
    Nmax = int(np.ceil((Ntot / Ncrit) + Ndisk * Nmin))
    RR = [0] * Nmax
    Lz = [0] * Nmax

    # memory allocation for energy
    TT = [0] * (Ndisk + 1)
    WW = [0] * (Ndisk + 1)

    # set partition for each component
    Ndat = [0] * Ndisk
    head = [0] * Ndisk
    head[0] = 0

    for diskID in range(Ndisk):
        kk = diskID + Nsphe

        # read attributes
        folder = "data" + str(kk) + "/"
        num = h5file[folder].attrs["number"][0]

        nunit = Ncrit
        if num < (Ncrit * Nmin):
            nunit = int(np.floor(num / Nmin))
        Ndat[diskID] = int(np.floor(num / nunit))

        # read particle position and velocity
        pos = h5file[folder + "position"]
        vel = h5file[folder + "velocity"]

        posx = pos[:, 0]
        posy = pos[:, 1]
        posz = pos[:, 2]
        velx = vel[:, 0]
        vely = vel[:, 1]
        velz = vel[:, 2]

        pot  = h5file[folder + "potential"]
        mass = h5file[folder + "mass"]

        WW[1 + diskID] = 0.5 * np.sum(mass * pot)

        # sort by R
        _R = posx * posx + posy * posy
        idx = np.argsort(_R)


        # analyze angular momentum profile of disk component(s)
        xx = np.empty(nunit)
        yy = np.empty(nunit)
        vx = np.empty(nunit)
        vy = np.empty(nunit)
        mm = np.empty(nunit)
        for ii in range(Ndat[diskID]):
            for jj in range(nunit):
                ll = idx[ii * nunit + jj]
                xx[jj] = posx[ll]
                yy[jj] = posy[ll]
                vx[jj] = velx[ll]
                vy[jj] = vely[ll]
                mm[jj] = mass[ll]

            R2 = xx * xx + yy * yy
            invR = 1.0 / np.sqrt(R2)

            # estimate z-component of total angular momentum Lz
            _RR = R2 * invR
            _lz = xx * vy - yy * vx
            # estimate mean rotation velocity v_rot
            vrot = np.mean(_lz / _RR)

            RR[head[diskID] + ii] = np.median(_RR)
            Lz[head[diskID] + ii] = np.median(_lz * mm)

            TT[1 + diskID] += 0.5 * np.sum(mm) * vrot * vrot


        if diskID != (Ndisk - 1):
            head[diskID + 1] = head[diskID] + Ndat[diskID]


        TT[0] += TT[1 + diskID]
        WW[0] += WW[1 + diskID]

    for kk in range(Nsphe):
        folder = "data" + str(kk) + "/"

        pot  = h5file[folder + "potential"]
        mass = h5file[folder + "mass"]

        WW[0] += 0.5 * np.sum(mass * pot)

    # close the HDF5 file
    h5file.close()

    tt = TT / np.abs(WW)


    for kk in range(Ndisk - 1, -1, -1):
        ax[0].plot(RR0[head0[kk] : head0[kk] + Ndat0[kk]], Lz0[head0[kk] : head0[kk] + Ndat0[kk]], linestyle = ls[kk], color = col[kk], label = tag[1 + kk] + " (" + r"$t = {:.3f}$ {:<}".format(0.0, "Gyr") + ")")



    # plot Lz profile of disk component(s)
    for kk in range(Ndisk - 1, -1, -1):
        ax[0].plot(RR[head[kk] : head[kk] + Ndat[kk]], Lz[head[kk] : head[kk] + Ndat[kk]], pt[kk], color = col[kk], label = tag[1 + kk] + " (" + r"$t = {:.3f}$ {:<}".format(time / 1000, "Gyr") + ")", markerfacecolor = "none")


    # set plot range
    ax[0].set_xlim([Rmin, Rmax])
    ax[0].set_ylim([Lmin, Lmax])
    # ax[0].semilogy()
    # ax[0].semilogx()
    ax[0].loglog()
    ax[0].grid()

    # set label
    ax[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode('UTF-8')))
    ax[0].set_ylabel(r"$L_z$ ($M_\odot$ {:<}".format(length_unit.decode('UTF-8')) + r" km s$^{-1}$)")

    # set legend
    handles, labels = ax[0].get_legend_handles_labels()
    # ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, prop = {'size' : 14}, loc = 'best')
    ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best')

    # output figure
    figname = "fig/" + filename + "_Lz_" + snapshot
    plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")


    # clear figure
    ax[0].cla()

    return time, tt



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


# read fundamental property of the system
target = "dat/" + filename + ".profile.h5"
data_file = h5py.File(target, "r")
mass_unit = data_file["/"].attrs["mass_astro_unit_name"][0]
length_unit = data_file["/"].attrs["length_astro_unit_name"][0]
velocity_unit = data_file["/"].attrs["velocity_astro_unit_name"][0]
time_unit = data_file["/"].attrs["time_astro_unit_name"][0]
Nkind = data_file["/"].attrs["kinds"][0]
data_file.close()


# read fundamental property of disk component(s)
target = "dat/" + filename + ".disk.h5"
if path.isfile(target) == True:
    disk_file = h5py.File(target, "r")
    Ndisk = disk_file["/"].attrs["kinds"][0]
    disk_file.close()
else:
    print("{:<} was not found.".format(target))
    sys.exit()


Nsphe = Nkind - Ndisk


Ndat0 = [0] * Ndisk
head0 = [0] * Ndisk
head0[0] = 0
h5file = h5py.File("dat/" + filename + ".split000.h5", "r")
Ntot = h5file["/"].attrs["number"][0]
Nmax = int(np.ceil((Ntot / Ncrit) + Ndisk * Nmin))
RR0 = [0] * Nmax
Lz0 = [0] * Nmax
for diskID in range(Ndisk):
    folder = "data" + str(diskID + Nsphe) + "/"

    num = h5file[folder].attrs["number"][0]
    nunit = Ncrit
    if num < (Ncrit * Nmin):
        nunit = int(np.floor(num / Nmin))
    Ndat0[diskID] = int(np.floor(num / nunit))

    # read particle position and velocity
    pos = h5file[folder + "position"]
    vel = h5file[folder + "velocity"]

    mass = h5file[folder + "mass"]
    posx = pos[:, 0]
    posy = pos[:, 1]
    velx = vel[:, 0]
    vely = vel[:, 1]

    _R = posx * posx + posy * posy
    idx = np.argsort(_R)

    xx = np.empty(nunit)
    yy = np.empty(nunit)
    vx = np.empty(nunit)
    vy = np.empty(nunit)
    mm = np.empty(nunit)
    for ii in range(Ndat0[diskID]):
        for jj in range(nunit):
            ll = idx[ii * nunit + jj]
            xx[jj] = posx[ll]
            yy[jj] = posy[ll]
            vx[jj] = velx[ll]
            vy[jj] = vely[ll]
            mm[jj] = mass[ll]

        R2 = xx * xx + yy * yy
        invR = 1.0 / np.sqrt(R2)

        _RR = R2 * invR
        _lz = xx * vy - yy * vx

        RR0[head0[diskID] + ii] = np.median(_RR)
        Lz0[head0[diskID] + ii] = np.median(_lz * mm)

    if diskID != (Ndisk - 1):
        head0[diskID + 1] = head0[diskID] + Ndat0[diskID]
h5file.close()




cores = mp.cpu_count()
# cores = int(np.ceil(cores / 2))

pool = mp.Pool(cores)
# ZZ = pool.map(make_map, range(nxpanel * nypanel))
args = [(ii, Nkind, Nsphe, Ndisk, mass_unit, length_unit, velocity_unit, time_unit, head0, Ndat0, RR0, Lz0) for ii in range(init, last + 1, 1)]
# time, tt = pool.map(wrapper, args)
ret = pool.map(wrapper, args)
pool.close()

Nfile = last + 1 - init
time = [0] * Nfile
tt = [[0 for ii in range(Nfile)] for jj in range(Ndisk + 1)]
for ii in range(Nfile):
    time[ii] = ret[ii][0]
    for jj in range(Ndisk + 1):
        tt[jj][ii] = ret[ii][1][jj]

fig = utils.set_figure(nxpanel, nypanel)
ax = [0] * nxpanel * nypanel
utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

for ii in range(Ndisk + 1 - 1, -1, -1):
    ax[0].plot(time, tt[ii], linestyle = ls[ii], marker = pt[ii], color = col[ii], label = tag[ii])

# set plot range
ax[0].set_xlim([time[0], time[Nfile - 1]])
ax[0].set_ylim([tmin, tmax])
# ax[0].semilogy()
# ax[0].semilogx()
# ax[0].loglog()
ax[0].grid()

# set label
ax[0].set_xlabel(r"$t$ ({:<})".format(time_unit.decode('UTF-8')))
ax[0].set_ylabel(r"$T / |W|$")

# set legend
handles, labels = ax[0].get_legend_handles_labels()
ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = "best")

# output figure
figname = "fig/" + filename + "_barmode"
plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")
