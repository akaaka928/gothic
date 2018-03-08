import numpy as np
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import os.path as path
import multiprocessing as mp

import utils as utils


series = "halocusp"
runname = series + "_run"
Ndark = 1
init = 0
# last = 0
last = 140
rmin, rmax = 1.0e-3, 5.0e+1
vmin, vmax = 0.0, 39.0
rhomin, rhomax = 2.0e+1, 1.0e+12
Ncrit = 2048 # number of particles for estimating physical quantities
# Ncrit = 8192 # number of particles for estimating physical quantities
Nmin = 16 # minimum data points for visualization
mass_astro2com = 1.0e-8
density_astro2com = 1.477569e+3

# set number of panels
nxpanel, nypanel = 1, 2

col = ["black", "red", "blue", "magenta"]


def draw_figures(Nkind, fileid, Nanal, radius, inirho, sigmar, simal, length_unit, velocity_unit, time_unit):
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    # for fileid in range(init, last + 1, 1):
    # get snapshot ID
    snapshot = "{:03d}".format(fileid)

    # get file name
    input_file = "dat/" + runname + ".split" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    Ntot = h5file["/"].attrs["number"][0]
    current_time = h5file["/"].attrs["time"][0]


    # memory allocation for velocity dispersion profile
    Nmax = int(np.ceil((Ntot / Ncrit) + Nkind * Nmin))
    rad = [0] * Nmax
    rho = [0] * Nmax
    sig = [0] * Nmax
    Rad = [0] * Nmax
    los = [0] * Nmax

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
            mass = h5file[folder + "mass"].value
            pos = h5file[folder + "position"].value
            vel = h5file[folder + "velocity"].value

            posx = pos[:, 0]
            posy = pos[:, 1]
            posz = pos[:, 2]
            velx = vel[:, 0]
            vely = vel[:, 1]
            velz = vel[:, 2]

            # find half-mass radius and center-of-mass of particles within the half-mass radius
            cx0, cy0, cz0 = 0.0, 0.0, 0.0
            vx0, vy0, vz0 = 0.0, 0.0, 0.0
            converge = False
            while converge == False:
                rrr = (posx - cx0) ** 2 + (posy - cy0) ** 2 + (posz - cz0) ** 2
                idx = np.argsort(rrr)

                mtot = 0.0
                comx, comy, comz = 0.0, 0.0, 0.0
                vx0 = 0.0
                vy0 = 0.0
                vz0 = 0.0
                for ii in range(int(num / 2)):
                    jj = idx[ii]
                    mm = mass[jj]
                    mtot += mm
                    comx += mm * posx[jj]
                    comy += mm * posy[jj]
                    comz += mm * posz[jj]
                    vx0  += mm * velx[jj]
                    vy0  += mm * vely[jj]
                    vz0  += mm * velz[jj]
                comx /= mtot
                comy /= mtot
                comz /= mtot
                vx0  /= mtot
                vy0  /= mtot
                vz0  /= mtot

                # convergence check
                dr2 = (cx0 - comx) ** 2 + (cy0 - comy) ** 2 + (cz0 - comz) ** 2
                cx0 = comx
                cy0 = comy
                cz0 = comz
                if dr2 < 1.0e-6:
                    converge = True

            # sort from radius from the center-of-mass
            hor = (posx - cx0) ** 2 + (posy - cy0) ** 2
            idx = np.argsort(hor)

            # analyze line-of-sight velocity dispersion profile of all component(s)
            xx = np.empty(nunit)
            yy = np.empty(nunit)
            zz = np.empty(nunit)
            vx = np.empty(nunit)
            vy = np.empty(nunit)
            vz = np.empty(nunit)
            for ii in range(Ndat[kk]):
                for jj in range(nunit):
                    mm = idx[ii * nunit + jj]
                    xx[jj] = posx[mm] - cx0
                    yy[jj] = posy[mm] - cy0
                    vz[jj] = velz[mm] - vz0

                R2 = xx * xx + yy * yy
                invR = 1.0 / np.sqrt(R2)

                Rad[head[kk] + ii] = np.median(R2 * invR)
                los[head[kk] + ii] = np.std(vz)


            # sort from radius from the center-of-mass
            rrr = (posx - cx0) ** 2 + (posy - cy0) ** 2 + (posz - cz0) ** 2
            idx = np.argsort(rrr)

            # analyze volume-density profile and velocity-dispersion profile
            for ii in range(Ndat[kk]):
                mtot = 0.0
                for jj in range(nunit):
                    mm = idx[ii * nunit + jj]
                    xx[jj] = posx[mm] - cx0
                    yy[jj] = posy[mm] - cy0
                    zz[jj] = posz[mm] - cz0
                    vx[jj] = velx[mm] - vx0
                    vy[jj] = vely[mm] - vy0
                    vz[jj] = velz[mm] - vz0
                    mtot += mass[mm]

                r2 = xx * xx + yy * yy + zz * zz
                invr = 1.0 / np.sqrt(r2)

                vr = (xx * vx + yy * vy + zz * vz) * invr

                rad[head[kk] + ii] = np.median(r2 * invr)
                rho[head[kk] + ii] = mtot / (4.0 * np.pi / 3.0 * (r2[nunit - 1] ** 1.5 - r2[0] ** 1.5))
                sig[head[kk] + ii] = np.std(vr)


            if kk != (Nkind - 1):
                head[kk + 1] = head[kk] + Ndat[kk]
        else:
            Ndat[kk] = 0
            if kk != (Nkind - 1):
                head[kk + 1] = head[kk] + Ndat[kk]

    # close the HDF5 file
    h5file.close()


    # plot velocity dispersion profile of all component(s) (analytic curve)
    for kk in range(Nanal - 1, -1, -1):
        if Ndat[kk] > 0:
            ax[0].plot(radius[kk], sigmal[kk], linestyle = ":", color = col[kk], label = r"$\sigma_\mathrm{los}$" + " (ID {:<})".format(kk))
            ax[0].plot(radius[kk], sigmar[kk], linestyle = "-", color = col[kk], label = r"$\sigma_r$" + " (ID {:<})".format(kk))

            ax[1].plot(radius[kk], inirho[kk], linestyle = "-", color = col[kk], label = r"$\rho$" + " (ID {:<})".format(kk))

    # plot velocity dispersion profile of all component(s) (N-body particle)
    for kk in range(Nkind - 1, -1, -1):
        if Ndat[kk] > 0:
            ax[0].plot(Rad[head[kk] : head[kk] + Ndat[kk]], los[head[kk] : head[kk] + Ndat[kk]], "s", color = col[kk], label = r"$\sigma_\mathrm{los}$" + " (ID {:<})".format(kk))
            ax[0].plot(rad[head[kk] : head[kk] + Ndat[kk]], sig[head[kk] : head[kk] + Ndat[kk]], "o", color = col[kk], label = r"$\sigma_r$" + " (ID {:<})".format(kk))

            ax[1].plot(rad[head[kk] : head[kk] + Ndat[kk]], rho[head[kk] : head[kk] + Ndat[kk]], "o", color = col[kk], label = r"$\rho$" + " (ID {:<})".format(kk))


    # set plot range
    ax[0].set_ylim([vmin, vmax])
    ax[0].semilogx()
    ax[1].set_ylim([rhomin, rhomax])
    ax[1].loglog()

    for ii in range(nypanel):
        ax[ii].set_xlim([rmin, rmax])
        ax[ii].grid()

    # set label
    ax[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode('UTF-8')))
    # ax[0].set_ylabel(r"$\sigma$ ({:<})".format(velocity_unit.decode('UTF-8')))
    ax[0].set_ylabel(r"$\sigma$ (km s$^{-1}$)")
    ax[1].set_ylabel(r"$\rho$ ($M_\odot$ {:<}".format(length_unit.decode('UTF-8')) + r"$^{-3}$)")

    # set legend
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best')

    # add current time
    fig.suptitle(r"$t = {:.3f}$ {:<}".format(current_time / 1000, "Gyr"))

    # output figure
    figname = "fig/" + runname + "_radial" + snapshot
    # plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".png", format = "png", dpi = 150, bbox_inches = "tight")

    # clear figure
    ax[0].cla()


def wrapper(argv):
    return draw_figures(*argv)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 28

rho2astro = density_astro2com / mass_astro2com

# read analytic profile of all component(s) in host component
target = "dat/" + series + ".profile.h5"
data_file = h5py.File(target, "r")
length_unit = data_file["/"].attrs["length_astro_unit_name"][0]
velocity_unit = data_file["/"].attrs["velocity_astro_unit_name"][0]
time_unit = data_file["/"].attrs["time_astro_unit_name"][0]
Nanal = data_file["/"].attrs["kinds"][0]
Ndata = data_file["/data0/"].attrs["num"][0]
radius = [0] * Nanal * Ndata
inirho = [0] * Nanal * Ndata
sigmar = [0] * Nanal * Ndata
sigmal = [0] * Nanal * Ndata
for kk in range(Nanal):
    folder = "data" + str(kk) + "/"
    radius[kk] = data_file[folder + "rad"].value
    inirho[kk] = data_file[folder + "rho"].value
    sigmar[kk] = data_file[folder + "sigma_r"].value
    sigmal[kk] = data_file[folder + "sigma_los"].value
    # unit conversion for density profile
    for ii in range(Ndata):
        inirho[kk][ii] *= rho2astro
data_file.close()


# read number of all component(s)
txtfile = open("doc/" + runname + ".summary.txt", "r")
unit = txtfile.readline()
tmp = txtfile.readline()
Nkind = int(tmp[0])
# Nsphe = int(tmp[2])
txtfile.close()


# cores = mp.cpu_count()
cores = int(np.ceil(mp.cpu_count() / 2))

pool = mp.Pool(cores)
args = [(Nkind, ii, Nanal, radius, inirho, sigmar, sigmal, length_unit, velocity_unit, time_unit) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
