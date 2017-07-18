import numpy as np
import h5py

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt

import os.path as path

import multiprocessing as mp

import utils as utils


# specify plot target
filename = "ltg"
# filename = "spherical"
init = 0
last = 47
# last = 0
Ncrit = 2048 # number of particles for estimating physical quantities
# Ncrit = 16384 # number of particles for estimating physical quantities
Nmin = 16 # minimum data points for visualization


# set plot range
rrmin, rrmax = 1.0e-1, 2.5e+2
vrmin, vrmax = 0.0, 250.0
RRmin, RRmax = 0.0, 50.0
vRmin, vRmax = 0.0, 70.0


# set number of panels
nxpanel, nypanel = 1, 1


col = ["black", "red", "blue", "magenta"]


def draw_figures(fileid, Nkind, radius, sigmar, sigmal, Nsphe, Ndisk, Radius, sigmaR, sigmap, sigmaz, length_unit, velocity_unit, time_unit):
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
    current_time = h5file["/"].attrs["time"][0]


    # memory allocation for velocity dispersion profile
    Nmax      = int(np.ceil((Ntot / Ncrit) + Nkind * Nmin))
    Nmax_disk = int(np.ceil((Ntot / Ncrit) + Ndisk * Nmin))
    rr = [0] * Nmax
    sr = [0] * Nmax
    RR = [0] * Nmax
    sz = [0] * Nmax
    sR = [0] * Nmax_disk
    sp = [0] * Nmax_disk

    # set partition for each component
    head = [0] * Nkind
    Ndat = [0] * Nkind
    head[0] = 0
    if Ndisk > 0:
        head_disk = [0] * Ndisk
        Ndat_disk = [0] * Ndisk
        head_disk[0] = 0

    for kk in range(Nkind):
        disk = False
        diskID = 0
        if kk >= Nsphe:
            disk = True
            diskID = kk - Nsphe


        # read attributes
        folder = "data" + str(kk) + "/"
        num = h5file[folder].attrs["number"][0]

        if num > 1:
            # skip analysis for central black hole
            nunit = Ncrit
            if num < (Ncrit * Nmin):
                nunit = int(np.floor(num / Nmin))

            Ndat[kk] = int(np.floor(num / nunit))
            if disk == True:
                Ndat_disk[diskID] = int(np.floor(num / nunit))

            # read particle position and velocity
            pos = h5file[folder + "position"].value
            vel = h5file[folder + "velocity"].value

            posx = pos[:, 0]
            posy = pos[:, 1]
            posz = pos[:, 2]
            velx = vel[:, 0]
            vely = vel[:, 1]
            velz = vel[:, 2]

            # sort by R
            _R = posx * posx + posy * posy
            idx = np.argsort(_R)

            # analyze velocity dispersion profile
            xx = np.empty(nunit)
            yy = np.empty(nunit)
            zz = np.empty(nunit)
            vx = np.empty(nunit)
            vy = np.empty(nunit)
            vz = np.empty(nunit)
            for ii in range(Ndat[kk]):
                for jj in range(nunit):
                    ll = idx[ii * nunit + jj]
                    xx[jj] = posx[ll]
                    yy[jj] = posy[ll]
                    zz[jj] = posz[ll]
                    vx[jj] = velx[ll]
                    vy[jj] = vely[ll]
                    vz[jj] = velz[ll]

                R2 = xx * xx + yy * yy
                invR = 1.0 / np.sqrt(R2)

                RR[head[kk] + ii] = np.median(R2 * invR)
                sz[head[kk] + ii] = np.std(vz)

                if disk == True:
                    vR = ( xx * vx + yy * vy) * invR
                    vp = (-yy * vx + xx * vz) * invR
                    sR[head_disk[diskID] + ii] = np.std(vR)
                    sp[head_disk[diskID] + ii] = np.std(vp)


            # sort by r
            _r = posx * posx + posy * posy + posz * posz
            idx = np.argsort(_r)

            # analyze velocity dispersion profile
            for ii in range(Ndat[kk]):
                for jj in range(nunit):
                    ll = idx[ii * nunit + jj]
                    xx[jj] = posx[ll]
                    yy[jj] = posy[ll]
                    zz[jj] = posz[ll]
                    vx[jj] = velx[ll]
                    vy[jj] = vely[ll]
                    vz[jj] = velz[ll]

                r2 = xx * xx + yy * yy + zz * zz
                invr = 1.0 / np.sqrt(r2)

                vr = (xx * vx + yy * vy + zz * vz) * invr

                rr[head[kk] + ii] = np.median(r2 * invr)
                sr[head[kk] + ii] = np.std(vr)


            if kk != (Nkind - 1):
                head[kk + 1] = head[kk] + Ndat[kk]
            if disk == True:
                if diskID != (Ndisk - 1):
                    head_disk[diskID + 1] = head_disk[diskID] + Ndat_disk[diskID]

        else:
            Ndat[kk] = 0
            if kk != (Nkind - 1):
                head[kk + 1] = head[kk] + Ndat[kk]

    # close the HDF5 file
    h5file.close()


    if Ndisk > 0:
        # plot velocity dispersion profile of disk component(s) (analytic curve)
        for kk in range(Ndisk - 1, -1, -1):
            # ax[0].plot(Radius[kk], sigmap[kk], linestyle = "--", color = col[kk], label = r"$\sigma_p$" + " (disk {:<})".format(kk))
            ax[0].plot(Radius[kk], sigmaz[kk], linestyle = ":", color = col[kk], label = r"$\sigma_z$" + " (disk {:<})".format(kk))
            ax[0].plot(Radius[kk], sigmaR[kk], linestyle = "-", color = col[kk], label = r"$\sigma_R$" + " (disk {:<})".format(kk))

        # plot velocity dispersion profile of disk component(s) (N-body particle)
        for kk in range(Ndisk - 1, -1, -1):
            # ax[0].plot(RR[head_disk[kk] : head_disk[kk] + Ndat_disk[kk]], sp[head_disk[kk] : head_disk[kk] + Ndat_disk[kk]], "D", color = col[kk], label = r"$\sigma_p$" + " (disk {:<})".format(kk))
            ax[0].plot(RR[head[Nsphe + kk] : head[Nsphe + kk] + Ndat[Nsphe + kk]], sz[head[Nsphe + kk] : head[Nsphe + kk] + Ndat[Nsphe + kk]], "s", color = col[kk], label = r"$\sigma_z$" + " (disk {:<})".format(kk))
            ax[0].plot(RR[head[Nsphe + kk] : head[Nsphe + kk] + Ndat[Nsphe + kk]], sR[head_disk[kk] : head_disk[kk] + Ndat_disk[kk]], "o", color = col[kk], label = r"$\sigma_R$" + " (disk {:<})".format(kk))


        # set plot range
        ax[0].set_xlim([RRmin, RRmax])
        ax[0].set_ylim([vRmin, vRmax])
        ax[0].grid()

        # set label
        ax[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode('UTF-8')))
        ax[0].set_ylabel(r"$\sigma$ ({:<})".format(velocity_unit.decode('UTF-8')))

        # set legend
        handles, labels = ax[0].get_legend_handles_labels()
        ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, prop = {'size' : 14}, loc = 'best')

        # output figure
        figname = "fig/" + filename + "_vdisp_disk_" + snapshot
        plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")


        # clear figure
        ax[0].cla()


    # plot velocity dispersion profile of all component(s) (analytic curve)
    for kk in range(Nsphe - 1, -1, -1):
        if Ndat[kk] > 0:
            ax[0].plot(radius[kk], sigmal[kk], linestyle = ":", color = col[kk], label = r"$\sigma_\mathrm{los}$" + " (component {:<})".format(kk))
            ax[0].plot(radius[kk], sigmar[kk], linestyle = "-", color = col[kk], label = r"$\sigma_r$" + " (component {:<})".format(kk))

    # plot velocity dispersion profile of all component(s) (N-body particle)
    for kk in range(Nsphe - 1, -1, -1):
        if Ndat[kk] > 0:
            ax[0].plot(RR[head[kk] : head[kk] + Ndat[kk]], sz[head[kk] : head[kk] + Ndat[kk]], "s", color = col[kk], label = r"$\sigma_\mathrm{los}$" + " (component {:<})".format(kk))
            ax[0].plot(rr[head[kk] : head[kk] + Ndat[kk]], sr[head[kk] : head[kk] + Ndat[kk]], "o", color = col[kk], label = r"$\sigma_r$" + " (component {:<})".format(kk))


    # set plot range
    ax[0].set_xlim([rrmin, rrmax])
    ax[0].set_ylim([vrmin, vrmax])
    ax[0].semilogx()
    ax[0].grid()

    # set label
    ax[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode('UTF-8')))
    ax[0].set_ylabel(r"$\sigma$ ({:<})".format(velocity_unit.decode('UTF-8')))

    # set legend
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, prop = {'size' : 14}, loc = 'best')

    # output figure
    figname = "fig/" + filename + "_vdisp_" + snapshot
    plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")

    # clear figure
    ax[0].cla()


def wrapper(argv):
    return draw_figures(*argv)






# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
plt.rcParams['font.size'] = 16


# read analytic profile of all component(s)
target = "dat/" + filename + ".profile.h5"
data_file = h5py.File(target, "r")
length_unit = data_file["/"].attrs["length_astro_unit_name"][0]
velocity_unit = data_file["/"].attrs["velocity_astro_unit_name"][0]
time_unit = data_file["/"].attrs["time_astro_unit_name"][0]
Nkind = data_file["/"].attrs["kinds"][0]
Ndata = data_file["/data0/"].attrs["num"][0]
radius = [0] * Nkind * Ndata
sigmar = [0] * Nkind * Ndata
sigmal = [0] * Nkind * Ndata
for kk in range(Nkind):
    folder = "data" + str(kk) + "/"
    radius[kk] = data_file[folder + "rad"].value
    sigmar[kk] = data_file[folder + "sigma_r"].value
    sigmal[kk] = data_file[folder + "sigma_los"].value
data_file.close()


# read analytic profile for disk component(s)
target = "dat/" + filename + ".disk.h5"
if path.isfile(target) == True:
    disk_file = h5py.File(target, "r")
    Ndisk = disk_file["/"].attrs["kinds"][0]
    Ndata = disk_file["/data0/1D_data/"].attrs["num"][0]
    Radius = [0] * Ndisk * Ndata
    sigmaR = [0] * Ndisk * Ndata
    sigmap = [0] * Ndisk * Ndata
    sigmaz = [0] * Ndisk * Ndata
    for kk in range(Ndisk):
        folder = "data" + str(kk) + "/1D_data/"
        Radius[kk] = disk_file[folder + "radius"].value
        sigmaR[kk] = disk_file[folder + "sigmaR"].value
        sigmap[kk] = disk_file[folder + "sigmap"].value
        sigmaz[kk] = disk_file[folder + "sigmaz"].value
    disk_file.close()
else:
    Ndisk = 0
    Radius = [0]
    sigmaR = [0]
    sigmap = [0]
    sigmaz = [0]


Nsphe = Nkind - Ndisk


# for fileid in range(init, last + 1, 1):
#     draw_figures(fileid, Nkind, radius, sigmar, sigmal, Nsphe, Ndisk, Radius, sigmaR, sigmap, sigmaz, length_unit, velocity_unit, time_unit)


cores = mp.cpu_count()
# cores = int(np.ceil(cores / 2))

pool = mp.Pool(cores)
# ZZ = pool.map(make_map, range(nxpanel * nypanel))
args = [(ii, Nkind, radius, sigmar, sigmal, Nsphe, Ndisk, Radius, sigmaR, sigmap, sigmaz, length_unit, velocity_unit, time_unit) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
