import numpy as np
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import multiprocessing as mp

import utils as utils


# specify plot target
dir_nbody = "/home/ymiki/calculation/nbody/tree/dat/"
dir_fixed = "dat/"
file_nbody = "ltg"
file_fixed = "ltg_disk"
init = 0
last = 3
hidx_disk_nbody = 8007308

# set plot range
rmin, rmax = 1.0e-2, 1.0e+2
Emin, Emax = 1.0e-3, 1.0e+2

# set number of panels
nxpanel, nypanel = 2, 2

col = ["black", "red", "blue", "magenta"]


def draw_figures(fileid):
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    # for fileid in range(init, last + 1, 1):
    # get snapshot ID
    snapshot = "{:03d}".format(fileid)

    # get file name
    input_nbody = dir_nbody + file_nbody + ".snp" + snapshot + ".h5"
    input_fixed = dir_fixed + file_fixed + ".snp" + snapshot + ".h5"

    # read snapshot
    h5file_nbody = h5py.File(input_nbody, "r")
    h5file_fixed = h5py.File(input_fixed, "r")

    # read attributes
    Nfull = h5file_nbody["/snp/"].attrs["number"][0]
    Npart = h5file_fixed["/snp/"].attrs["number"][0]
    length_unit = h5file_nbody["/unit_system/axis labels"].attrs["length_astro_unit_name"][0]

    if (Npart + hidx_disk_nbody) != Nfull:
        alert

    # memory allocation for error analysis
    idx = [0] * Npart
    rad = [0] * Npart
    err_Phi = [0] * Npart
    err_ar = [0] * Npart
    err_az = [0] * Npart
    err_aR = [0] * Npart

    # memory allocation for comparing radial profiles
    rad_nbody = [0] * Nfull
    rad_fixed = [0] * Npart
    Phi_nbody = [0] * Nfull
    Phi_local = [0] * Npart
    Phi_field = [0] * Npart
    acc_nbody = [0] * Nfull
    acc_local = [0] * Npart
    acc_field = [0] * Npart

    # read particle position and velocity
    folder = "snp/"
    pos_part = h5file_fixed[folder + "position"]
    acc_part = h5file_fixed[folder + "acceleration"]
    pot_part = h5file_fixed[folder + "potential"]
    idx_part = h5file_fixed[folder + "index"]
    acc_ext_part = h5file_fixed[folder + "acceleration_external"]
    pot_ext_part = h5file_fixed[folder + "potential_external"]
    pos_full = h5file_nbody[folder + "position"]
    acc_full = h5file_nbody[folder + "acceleration"]
    pot_full = h5file_nbody[folder + "potential"]
    idx_full = h5file_nbody[folder + "index"]

    ax_full = acc_full[:, 0]
    ay_full = acc_full[:, 1]
    az_full = acc_full[:, 2]

    ax_part = acc_part[:, 0]# + acc_ext_part[:, 0]
    ay_part = acc_part[:, 1]# + acc_ext_part[:, 1]
    az_part = acc_part[:, 2]# + acc_ext_part[:, 2]

    ax_ext_part = acc_ext_part[:, 0]
    ay_ext_part = acc_ext_part[:, 1]
    az_ext_part = acc_ext_part[:, 2]

    xx_full = pos_full[:, 0]
    yy_full = pos_full[:, 1]
    zz_full = pos_full[:, 2]

    xx_part = pos_part[:, 0]
    yy_part = pos_part[:, 1]
    zz_part = pos_part[:, 2]

    # close the HDF5 file
    h5file_nbody.close()
    h5file_fixed.close()

    # start error analysis
    src_part = np.argsort(idx_part)
    src_full = np.argsort(idx_full)

    for ii in range(Npart):
        idx_nbody = src_full[hidx_disk_nbody + ii]
        idx_fixed = src_part[                  ii]

        pot_nbody = pot_full[idx_nbody]
        pot_fixed = pot_part[idx_fixed] + pot_ext_part[idx_fixed]

        ax_nbody = ax_full[idx_nbody]
        ay_nbody = ay_full[idx_nbody]
        az_nbody = az_full[idx_nbody]

        ax_fixed = ax_part[idx_fixed] + ax_ext_part[idx_fixed]
        ay_fixed = ay_part[idx_fixed] + ay_ext_part[idx_fixed]
        az_fixed = az_part[idx_fixed] + az_ext_part[idx_fixed]

        xx_nbody = xx_full[idx_nbody]
        yy_nbody = yy_full[idx_nbody]
        zz_nbody = zz_full[idx_nbody]

        xx_fixed = xx_part[idx_fixed]
        yy_fixed = yy_part[idx_fixed]
        zz_fixed = zz_part[idx_fixed]

        idx[ii] = idx_fixed
        rad[ii] = np.sqrt(xx_fixed * xx_fixed + yy_fixed * yy_fixed + zz_fixed * zz_fixed)

        err_Phi[ii] = (pot_fixed - pot_nbody) / pot_nbody

        aR_nbody = (ax_nbody * xx_nbody + ay_nbody * yy_nbody) / np.sqrt(xx_nbody * xx_nbody + yy_nbody * yy_nbody)
        ar_nbody = (ax_nbody * xx_nbody + ay_nbody * yy_nbody + az_nbody * zz_nbody) / np.sqrt(xx_nbody * xx_nbody + yy_nbody * yy_nbody + zz_nbody * zz_nbody)

        aR_fixed = (ax_fixed * xx_fixed + ay_fixed * yy_fixed) / np.sqrt(xx_fixed * xx_fixed + yy_fixed * yy_fixed)
        ar_fixed = (ax_fixed * xx_fixed + ay_fixed * yy_fixed + az_fixed * zz_fixed) / np.sqrt(xx_fixed * xx_fixed + yy_fixed * yy_fixed + zz_fixed * zz_fixed)

        err_az[ii] = (az_fixed - az_nbody) / az_nbody
        err_ar[ii] = (ar_fixed - ar_nbody) / ar_nbody
        err_aR[ii] = (aR_fixed - aR_nbody) / aR_nbody


    # generate radial profiles
    rad_fixed = np.sqrt(xx_part * xx_part + yy_part * yy_part + zz_part * zz_part)
    rad_nbody = np.sqrt(xx_full * xx_full + yy_full * yy_full + zz_full * zz_full)

    src_nbody = np.argsort(rad_nbody)
    src_fixed = np.argsort(rad_fixed)

    rad_nbody = np.sort(rad_nbody)
    rad_fixed = np.sort(rad_fixed)

    for ii in range(Npart):
        idx_fixed = src_fixed[ii]

        Phi_local[ii] = pot_part[idx_fixed]
        Phi_field[ii] = pot_ext_part[idx_fixed]

        ax_local = ax_part[idx_fixed]
        ay_local = ay_part[idx_fixed]
        az_local = az_part[idx_fixed]
        ax_field = ax_ext_part[idx_fixed]
        ay_field = ay_ext_part[idx_fixed]
        az_field = az_ext_part[idx_fixed]

        acc_local[ii] = np.sqrt(ax_local * ax_local + ay_local * ay_local + az_local * az_local)
        acc_field[ii] = np.sqrt(ax_field * ax_field + ay_field * ay_field + az_field * az_field)

    for ii in range(Nfull):
        idx_nbody = src_nbody[ii]

        Phi_nbody[ii] = pot_full[idx_nbody]

        ax_nbody = ax_full[idx_nbody]
        ay_nbody = ay_full[idx_nbody]
        az_nbody = az_full[idx_nbody]

        acc_nbody[ii] = np.sqrt(ax_nbody * ax_nbody + ay_nbody * ay_nbody + az_nbody * az_nbody)


    # plot relative error in linear scale
    # plot relative error of potential as function of particle index
    ax[0].plot(idx, err_Phi, ",", color = "black", rasterized = True)
    ax[0].set_xlabel("Particle ID")
    ax[0].set_ylabel("Error of " + r"$\Phi(R, z)$")

    # plot relative error of acceleration as function of particle index
    ax[1].plot(idx, err_az, ",", color = "blue", rasterized = True, label = r"$a_z$")
    ax[1].plot(idx, err_aR, ",", color = "red", rasterized = True, label = r"$a_R$")
    ax[1].plot(idx, err_ar, ",", color = "black", rasterized = True, label = r"$a_r$")
    ax[1].set_ylabel("Error of " + r"$\bm{a}(R, z)$")

    # plot relative error of potential as function of particle position
    ax[2].plot(rad, err_Phi, ",", color = "black", rasterized = True)
    ax[2].set_xlim([rmin, rmax])
    ax[2].semilogx()
    ax[2].set_xlabel(r"$r$ ({:<})".format(length_unit.decode('UTF-8')))

    # plot relative error of acceleration as function of particle position
    ax[3].plot(rad, err_az, ",", color = "blue", rasterized = True)
    ax[3].plot(rad, err_aR, ",", color = "red", rasterized = True)
    ax[3].plot(rad, err_ar, ",", color = "black", rasterized = True)
    ax[3].set_xlim([rmin, rmax])
    ax[3].semilogx()

    for ii in range(nxpanel * nypanel):
        ax[ii].grid()
        ax[ii].set_ylim([-Emax, Emax])

    # # set legend
    # handles, labels = ax[1].get_legend_handles_labels()
    # ax[1].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best')

    # output figure
    figname = "fig/" + file_nbody + "_err_lin_" + snapshot
    plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")

    # # clear figure
    # for ii in range(nxpanel * nypanel):
    #     ax[ii].cla()

    # plot relative error in logarithmic scale
    for ii in range(nxpanel * nypanel):
        ax[ii].set_ylim([Emin, Emax])

    ax[0].semilogy()
    ax[1].semilogy()
    ax[2].loglog()
    ax[3].loglog()

    ax[0].text(0, Emin * 10, "med" + r"$(\Phi) = {:.1e}$".format(np.median(err_Phi)))
    ax[0].text(0, Emin * 2, "med" + r"$(|\Phi|) = {:.1e}$".format(np.median(np.abs(err_Phi))))

    ax[1].text(0, Emin * 50, "med" + r"$(|a_r|) = {:.1e}$".format(np.median(np.abs(err_ar))))
    ax[1].text(0, Emin * 10, "med" + r"$(|a_R|) = {:.1e}$".format(np.median(np.abs(err_aR))))
    ax[1].text(0, Emin * 2, "med" + r"$(|a_z|) = {:.1e}$".format(np.median(np.abs(err_az))))

    ax[3].text(rmin * 1.5, Emin * 50, "med" + r"$(a_r) = {:.1e}$".format(np.median(err_ar)))
    ax[3].text(rmin * 1.5, Emin * 10, "med" + r"$(a_R) = {:.1e}$".format(np.median(err_aR)))
    ax[3].text(rmin * 1.5, Emin * 2, "med" + r"$(a_z) = {:.1e}$".format(np.median(err_az)))


    # output figure
    figname = "fig/" + file_nbody + "_err_log_" + snapshot
    plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")

    # clear figure
    for ii in range(nxpanel * nypanel):
        ax[ii].cla()

    # plot radial profiles
    ax[0].plot(rad_nbody, Phi_nbody, "-", color = "black", rasterized = True, label = "full N-body")
    ax[0].plot(rad_fixed, Phi_local, "-", color = "red", rasterized = True, label = "N-body part")
    ax[0].plot(rad_fixed, Phi_field, "-", color = "blue", rasterized = True, label = "potential part")
    ax[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode('UTF-8')))
    ax[0].set_ylabel(r"$\Phi(R, z)$")
    ax[0].set_xlim([rmin, rmax])
    ax[0].semilogx()

    ax[1].plot(rad_nbody, acc_nbody, "-", color = "black", rasterized = True, label = "full N-body")
    ax[1].plot(rad_fixed, acc_local, "-", color = "red", rasterized = True, label = "N-body part")
    ax[1].plot(rad_fixed, acc_field, "-", color = "blue", rasterized = True, label = "potential part")
    ax[1].set_ylabel(r"$|a|$")
    ax[1].set_xlim([rmin, rmax])
    ax[1].loglog()

    # set legend
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best')

    # output figure
    figname = "fig/" + file_nbody + "_compare_" + snapshot
    plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")

    # clear figure
    for ii in range(nxpanel * nypanel):
        ax[ii].cla()





# def wrapper(argv):
#     return draw_figures(*argv)




# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
plt.rcParams['font.size'] = 28

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.rcParams['text.latex.preamble'] = [r"\usepackage{bm}"]

cores = mp.cpu_count()
pool = mp.Pool(cores)
# args = [(ii) for ii in range(init, last + 1, 1)]
# pool.map(wrapper, args)
pool.map(draw_figures, range(init, last + 1, 1))
pool.close()
