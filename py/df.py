import numpy as np
import h5py

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt

import os.path as path

import multiprocessing as mp

import utils as utils


# specify plot target
# filename = "ltg"
# filename = "spherical"
filename = "plummer"
init = 0
# last = 47
last = 0
# Ncrit = 128 # number of particles for estimating physical quantities
Ncrit = 2048 # number of particles for estimating physical quantities
# Ncrit = 8192 # number of particles for estimating physical quantities
# Ncrit = 16384 # number of particles for estimating physical quantities
Nmin = 16 # minimum data points for visualization


# set plot range
Emin, Emax = 0.0, 1.1
# Emin, Emax = 1.0e+4, 2.0e+6
# Emin, Emax = 1.0e+1, 2.0e+6
rmin, rmax = 1.0e-1, 2.5e+2
# fmin, fmax = 1.0e-14, 1.0e+18
fmin, fmax = 1.0e-6, 1.0e+0
vmin, vmax = -300.0, 300.0


# set number of panels
nxpanel, nypanel = 1, 1


col = ["black", "red", "blue", "magenta"]


def draw_figures(fileid, Nkind, rr, DF):
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    # for fileid in range(init, last + 1, 1):
    # get snapshot ID
    snapshot = "{:03d}".format(fileid)


    # generate Phi(r)
    input_file = "dat/" + filename + ".snp" + snapshot + ".h5"
    h5file = h5py.File(input_file, "r")

    # read attributes
    Ntot = h5file["/snp/"].attrs["number"][0]
    current_time = h5file["/snp/"].attrs["time"][0]

    potential = h5file["/snp/potential"].value
    position = h5file["/snp/position"].value
    radius = np.sqrt(position[:, 0] * position[:, 0] + position[:, 1] * position[:, 1] + position[:, 2] * position[:, 2])
    idx = np.argsort(radius)

    nunit = Ncrit
    if Ntot < (Ncrit * Nmin):
        nunit = int(np.floor(Ntot / Nmin))

    Ntbl = int(np.floor(Ntot / nunit))
    rad = [0] * Ntbl
    Phi = [0] * Ntbl

    kk = 0
    rad_sub = np.empty(nunit)
    pot_sub = np.empty(nunit)
    for ii in range(Ntbl):
        for jj in range(nunit):
            ll = idx[ii * nunit + jj]
            rad_sub[jj] = radius[ll]
            pot_sub[jj] = potential[ll]
        rad[kk] = np.median(rad_sub)
        Phi[kk] = np.median(pot_sub)
        kk += 1

    # close the HDF5 file
    h5file.close()

    # print(rad)
    # print(Phi)

    # get file name
    input_file = "dat/" + filename + ".split" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")




    # 1. generate Phi(r)
    # 2. then, get r_out


    # memory allocation for distribution function
    Nmax      = int(np.ceil((Ntot / Ncrit) + Nkind * Nmin))
    rout = [0] * Nmax
    val = [0] * Nmax

    # set partition for each component
    head = [0] * Nkind
    Ndat = [0] * Nkind
    head[0] = 0

    for kk in range(Nkind):
        # read attributes
        folder = "data" + str(kk) + "/"
        num = h5file[folder].attrs["number"][0]

        if num > 1:
            # skip analysis for central black hole
            nunit = Ncrit
            if num < (Ncrit * Nmin):
                nunit = int(np.floor(num / Nmin))

            Ndat[kk] = int(np.floor(num / nunit))

            # read particle position and velocity
            vel = h5file[folder + "velocity"].value
            pot = h5file[folder + "potential"].value
            mass = h5file[folder + "mass"].value

            # sort by apocenter
            v2 = vel[:, 0] * vel[:, 0] + vel[:, 1] * vel[:, 1] + vel[:, 2] * vel[:, 2]
            ene = pot + 0.5 * v2
            # print(ene)
            apo = [0] * num

            # find rad corresponding to ene == Phi using bisection
            for ii in range(num):
                il = 0
                ir = Ntbl - 1

                while( True ):
                    ic = int((il + ir) / 2)

                    if ((Phi[ic] - ene[ii]) * (Phi[il] - ene[ii]) <= 0.0):
                        ir = ic
                    else:
                        il = ic

                    if (1 + il) == ir:
                        break

                apo[ii] = rad[il] + (rad[ir] - rad[il]) * (ene[ii] - Phi[il]) / (Phi[ir] - Phi[il])

            idx = np.argsort(apo)

            # ax[0].hist(apo, bins = 100, normed = True)
            # fig.show()

            # extract distribution function
            msub = np.empty(nunit)
            rsub = np.empty(nunit)
            for ii in range(Ndat[kk]):
                for jj in range(nunit):
                    ll = idx[ii * nunit + jj]
                    msub[jj] = mass[ll]
                    rsub[jj] = apo[ll]
                rout[head[kk] + ii] = np.median(rsub)
                val[head[kk] + ii] = np.sum(msub) / (4.0 * np.pi * (rsub[nunit - 1] ** 3 - rsub[0] ** 3) / 3.0)


            # pos = h5file[folder + "position"].value
            # rad = np.sqrt(pos[:, 0] * pos[:, 0] + pos[:, 1] * pos[:, 1] + pos[:, 2] * pos[:, 2])
            # idx = np.argsort(rad)
            # msub = np.empty(nunit)
            # rsub = np.empty(nunit)
            # for ii in range(Ndat[kk]):
            #     for jj in range(nunit):
            #         ll = idx[ii * nunit + jj]
            #         rsub[jj] = rad[ll]
            #         msub[jj] = -pot[ll]
            #     rout[head[kk] + ii] = np.median(rsub)
            #     val[head[kk] + ii] = np.median(msub)

            # print(rout)
            # print(val)

            if kk != (Nkind - 1):
                head[kk + 1] = head[kk] + Ndat[kk]

        else:
            Ndat[kk] = 0
            if kk != (Nkind - 1):
                head[kk + 1] = head[kk] + Ndat[kk]

    # close the HDF5 file
    h5file.close()


    # plot distribution function of all component(s) (analytic curve)
    for kk in range(Nkind - 1, -1, -1):
        if Ndat[kk] > 0:
            ax[0].plot(rr[kk], DF[kk], linestyle = "-", color = col[kk], label = "analytic (component {:<})".format(kk))

    # plot distribution function of all component(s) (N-body particle)
    for kk in range(Nkind - 1, -1, -1):
        if Ndat[kk] > 0:
            ax[0].plot(rout[head[kk] : head[kk] + Ndat[kk]], val[head[kk] : head[kk] + Ndat[kk]], "o", color = col[kk], label = "particle (component {:<})".format(kk))


    # set plot range
    # ax[0].set_xlim([Emin, Emax])
    ax[0].set_xlim([rmin, rmax])
    ax[0].set_ylim([fmin, fmax])
    # ax[0].semilogy()
    ax[0].loglog()
    ax[0].grid()

    # # set label
    # ax[0].set_xlabel(r"$E$ ({:<})".format(length_unit.decode('UTF-8')))
    # ax[0].set_ylabel(r"$f(E)$ ({:<})".format(velocity_unit.decode('UTF-8')))

    # set legend
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, prop = {'size' : 14}, loc = 'best')

    # output figure
    figname = "fig/" + filename + "_df_" + snapshot
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


# we need transformation from sE to rad via potential
target = "dat/" + filename + ".profile.h5"
data_file = h5py.File(target, "r")
Nkind = data_file["/"].attrs["kinds"][0]
rad = data_file["/data0/rad"].value
psi = data_file["/data0/Psi"].value
for kk in range(1, Nkind):
    psi += data_file["/data" + str(kk) + "/Psi"].value

# read analytic distribution function
target = "dat/" + filename + ".df.h5"
data_file = h5py.File(target, "r")
# senergy_unit = data_file["/"].attrs["senergy_astro_unit_name"][0]
Ndata = data_file["/series0/"].attrs["num"][0]
sE = np.zeros((Nkind, Ndata))
DF = np.zeros((Nkind, Ndata))
rr = np.zeros((Nkind, Ndata))
for kk in range(Nkind):
    folder = "series" + str(kk) + "/"
    sE[kk] = data_file[folder + "energy"].value
    DF[kk] = data_file[folder + "DF"].value
    for ii in range(Ndata - 1):
        il = 0
        ir = rad.size - 1

        while( True ):
            ic = int((il + ir) / 2)

            if ((psi[ic] - sE[kk][ii]) * (psi[il] - sE[kk][ii]) <= 0.0):
                ir = ic
            else:
                il = ic

            if (1 + il) == ir:
                break

        rr[kk][ii] = rad[il] + (rad[ir] - rad[il]) * (sE[kk][ii] - psi[il]) / (psi[ir] - psi[il])

    rr[kk][Ndata - 1] = rr[kk][Ndata - 2]

data_file.close()

for fileid in range(init, last + 1, 1):
    draw_figures(fileid, Nkind, rr, DF)


# cores = mp.cpu_count()
# # cores = int(np.ceil(cores / 2))

# pool = mp.Pool(cores)
# args = [(ii, Nkind, sE, DF) for ii in range(init, last + 1, 1)]
# pool.map(wrapper, args)
# pool.close()
