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


osipkov_merritt = False


# filename = "mw"
# init = 0
# last = 127
# # last = 0
# lab = ["DM halo", "stellar halo", "bulge", "central BH", "thick disk", "thin disk"]

# filename = "m31"
# init = 0
# last = 47
# lab = ["DM halo", "stellar halo", "bulge", "disk"]

# filename = "k17disk"
# init = 0
# last = 399
# lab = ["halo", "bulge", "MBH", "disk"]

# # filename = "hd"
# filename = "disk"
# init = 0
# last = 13
# lab = ["halo", "disk"]

# # filename = "halocusp_run"
# # filename = "halocore1_run"
# filename = "halocore2_run"
# # filename = "halocore3_run"
# init = 0
# last = 140
# # last = 0
# lab = ["halo", "core", "GC"]

# filename = "hernquist"
# init = 0
# last = 47
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "hb"
# init = 0
# last = 47
# # last = 0
# lab = ["halo", "bulge"]
# osipkov_merritt = True

# filename = "etg"
# # Nskip = 0
# init = 0
# last = 47
# last = 15
# # last = 0
# lab = ["DM halo", "bulge", "stellar halo", "central BH"]

# filename = "satellite"
# init = 0
# last = 140
# # last = 0
# lab = ["bulge"]

# filename = "m12iso"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]

# filename = "m12ra2"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

filename = "m12ra1"
init = 0
last = 140
# last = 0
lab = ["halo"]
osipkov_merritt = True

# filename = "m12ra1_2"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "m12ra1_4"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "m12ra1_8"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "m09iso"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]

# filename = "m09ra2"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "m09ra1"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "m09ra1_2"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "m09ra1_4"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "m09ra1_8"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "a00iso"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]

# filename = "a00ra2"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "a00ra1"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "a00ra1_2"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True

# filename = "a00ra1_4"
# init = 0
# last = 140
# # last = 0
# lab = ["halo"]
# osipkov_merritt = True


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", ":", "-.", "--"]
Ncol = 8
col = ["black", "red", "blue", "magenta", "green", "brown", "cyan", "yellow"]

thicknessMeasure = 1.82


def draw_figure(fileid, kind, sphe, radmin, radmax, rhomin, rhomax, encmin, encmax, sigmin, sigmax, betmin, betmax, hormin, hormax, Sigmin, Sigmax, sigRmin, sigRmax, sigpmin, sigpmax, sigzmin, sigzmax, diskmin, diskmax, add_ini_sphe, sphe_rad, sphe_rho, sphe_enc, sphe_Sig, sphe_sig, sphe_sgt, sphe_bet):
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

    fig_rho = utils.set_figure(1, 1)
    fig_enc = utils.set_figure(1, 1)
    fig_sig = utils.set_figure(1, 1)
    fig_bet = utils.set_figure(1, 1)
    fig_Sig = utils.set_figure(1, 1)
    fig_ver = utils.set_figure(1, 1)
    fig_sig2D = utils.set_figure(1, 1)
    ax_rho = [0]
    ax_enc = [0]
    ax_sig = [0]
    ax_bet = [0]
    ax_Sig = [0]
    ax_ver = [0]
    ax_sig2D = [0]
    utils.locate_panels(fig_rho, ax_rho, 1, 1, True, True)
    utils.locate_panels(fig_enc, ax_enc, 1, 1, True, True)
    utils.locate_panels(fig_sig, ax_sig, 1, 1, True, True)
    utils.locate_panels(fig_bet, ax_bet, 1, 1, True, True)
    utils.locate_panels(fig_Sig, ax_Sig, 1, 1, True, True)
    utils.locate_panels(fig_ver, ax_ver, 1, 1, True, True)
    utils.locate_panels(fig_sig2D, ax_sig2D, 1, 1, True, True)


    for ii in range(kind):
        num = h5file["attr" + str(ii)].attrs["number"][0]
        if (num > 1):
            folder = "rad" + str(ii) + "/"
            rad = h5file[folder + "rad"]
            rho = h5file[folder + "rho"]
            enc = h5file[folder + "enc"]
            sig = h5file[folder + "sig"]
            if osipkov_merritt:
                sgt = h5file[folder + "tan"]
                bet = h5file[folder + "bet"]

            folder = "hor" + str(ii) + "/"
            hor = h5file[folder + "hor"]
            ver = h5file[folder + "height"]
            Sig = h5file[folder + "Sigma"]
            sigR = h5file[folder + "sigR"]
            sigp = h5file[folder + "sigp"]
            sigz = h5file[folder + "sigz"]

            # measure of the thickness, according to Rodionov & Sotnikova (2013)
            ver *= thicknessMeasure

            if add_ini_sphe:
                ax_rho[0].plot(sphe_rad[ii], sphe_rho[ii], linestyle = ls[ii % Nls], color = col[ii % Ncol])
                ax_enc[0].plot(sphe_rad[ii], sphe_enc[ii], linestyle = ls[ii % Nls], color = col[ii % Ncol])
                ax_Sig[0].plot(sphe_rad[ii], sphe_Sig[ii], linestyle = ls[ii % Nls], color = col[ii % Ncol])
                if osipkov_merritt:
                    ax_sig[0].plot(sphe_rad[ii], sphe_sig[ii], linestyle = ls[(2 * ii    ) % Nls], color = col[ii % Ncol], label = r"$\sigma_r$" + " (" + lab[ii] + ")")
                    ax_sig[0].plot(sphe_rad[ii], sphe_sgt[ii], linestyle = ls[(2 * ii + 1) % Nls], color = col[ii % Ncol], label = r"$\sigma_t$" + " (" + lab[ii] + ")")
                    ax_bet[0].plot(sphe_rad[ii], sphe_bet[ii], linestyle = ls[ii % Nls], color = col[ii % Ncol])
                else:
                    ax_sig[0].plot(sphe_rad[ii], sphe_sig[ii], linestyle = ls[ii % Nls], color = col[ii % Ncol])


            ax_rho[0].plot(rad, rho, pt[ii % Npt], linestyle = "None" if add_ini_sphe else ls[ii % Nls], color = col[ii % Ncol], label = lab[ii])
            ax_enc[0].plot(rad, enc, pt[ii % Npt], linestyle = "None" if add_ini_sphe else ls[ii % Nls], color = col[ii % Ncol], label = lab[ii])
            if osipkov_merritt:
                ax_bet[0].plot(rad, bet, pt[ii % Npt], linestyle = "None" if add_ini_sphe else ls[ii % Nls], color = col[ii % Ncol], label = lab[ii])
                ax_sig[0].plot(rad, sig, pt[ 0 % Npt], linestyle = "None" if add_ini_sphe else ls[ii % Nls], color = col[ii % Ncol], label = r"$\sigma_r$" + " (" + lab[ii] + ")")
                ax_sig[0].plot(rad, sgt, pt[ 1 % Npt], linestyle = "None" if add_ini_sphe else ls[ii % Nls], color = col[ii % Ncol], label = r"$\sigma_t$" + " (" + lab[ii] + ")")
            else:
                ax_sig[0].plot(rad, sig, pt[ii % Npt], linestyle = "None" if add_ini_sphe else ls[ii % Nls], color = col[ii % Ncol], label = lab[ii])
            ax_Sig[0].plot(hor, Sig, pt[ii % Npt], linestyle = "None" if add_ini_sphe else ls[ii % Nls], color = col[ii % Ncol], label = lab[ii])
            if ii >= sphe:
                ax_ver[0].plot(hor, ver, pt[(ii - sphe) % Npt], linestyle = ls[(ii - sphe) % Nls], color = col[(ii - sphe) % Ncol], label = lab[ii])
                ax_sig2D[0].plot(hor, sigR, pt[0 % Npt], linestyle = ls[(ii - sphe) % Nls], color = col[(ii - sphe) % Ncol], label = r"$\sigma_R$" + " (" + lab[ii] + ")")
                ax_sig2D[0].plot(hor, sigp, pt[1 % Npt], linestyle = ls[(ii - sphe) % Nls], color = col[(ii - sphe) % Ncol], label = r"$\sigma_p$" + " (" + lab[ii] + ")")
                ax_sig2D[0].plot(hor, sigz, pt[2 % Npt], linestyle = ls[(ii - sphe) % Nls], color = col[(ii - sphe) % Ncol], label = r"$\sigma_z$" + " (" + lab[ii] + ")")


    # close the HDF5 file
    h5file.close()


    ax_rho[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_enc[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_sig[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_Sig[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_ver[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode("UTF-8")))
    ax_sig2D[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode("UTF-8")))

    ax_rho[0].set_ylabel(r"$\rho$ ({:<})".format(density_unit.decode("UTF-8")))
    ax_enc[0].set_ylabel(r"$M_\mathrm{enc}$" + r" ({:<})".format(mass_unit.decode("UTF-8")))
    if osipkov_merritt:
        ax_sig[0].set_ylabel(r"$\sigma$ ({:<})".format(velocity_unit.decode("UTF-8")))
        ax_bet[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode("UTF-8")))
        ax_bet[0].set_ylabel(r"$\beta$")
    else:
        ax_sig[0].set_ylabel(r"$\sigma_r$ ({:<})".format(velocity_unit.decode("UTF-8")))
    ax_Sig[0].set_ylabel(r"$\Sigma$ ({:<})".format(col_density_unit.decode("UTF-8")))
    ax_ver[0].set_ylabel(r"${:.2f}$".format(thicknessMeasure) + r" $\mathrm{median}(|z|)$" + r" ({:<})".format(length_unit.decode("UTF-8")))
    ax_sig2D[0].set_ylabel(r"$\sigma$ ({:<})".format(velocity_unit.decode("UTF-8")))

    ax_rho[0].set_xlim(utils.scale_axis(radmin, radmax,  True))
    ax_rho[0].set_ylim(utils.scale_axis(rhomin, rhomax,  True))
    ax_enc[0].set_xlim(utils.scale_axis(radmin, radmax,  True))
    ax_enc[0].set_ylim(utils.scale_axis(encmin, encmax,  True))
    ax_sig[0].set_xlim(utils.scale_axis(radmin, radmax,  True))
    ax_sig[0].set_ylim(utils.scale_axis(sigmin, sigmax, False))

    ax_Sig[0].set_xlim(utils.scale_axis( hormin,  hormax, True))
    ax_Sig[0].set_ylim(utils.scale_axis( Sigmin,  Sigmax, True))
    ax_ver[0].set_xlim(utils.scale_axis( hormin,  hormax, True))
    ax_ver[0].set_ylim(utils.scale_axis(diskmin, thicknessMeasure * diskmax, False))
    ax_sig2D[0].set_xlim(utils.scale_axis( hormin,  hormax, True))
    ax_sig2D[0].set_ylim(utils.scale_axis(min([sigRmin, sigpmin, sigzmin]), max([sigRmax, sigpmax, sigzmax]), False))

    if osipkov_merritt:
        ax_bet[0].set_xlim(utils.scale_axis(radmin, radmax,  True))
        ax_bet[0].set_ylim(utils.scale_axis(betmin, betmax, False))

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

    if osipkov_merritt:
        ax_bet[0].semilogx()
        ax_bet[0].grid()

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

    if osipkov_merritt:
        handles, labels = ax_bet[0].get_legend_handles_labels()
        ax_bet[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')


    # add current time
    # fig_rho.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    # fig_enc.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    # fig_sig.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    # fig_Sig.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    # fig_ver.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    # fig_sig2D.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    # if osipkov_merritt:
    #     fig_bet.suptitle(r"$t = {:.3f}$ {:<}".format(time, time_unit.decode("UTF-8")))
    fig_rho.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})")
    fig_enc.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})")
    fig_sig.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})")
    fig_Sig.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})")
    fig_ver.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})")
    fig_sig2D.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})")
    if osipkov_merritt:
        fig_bet.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})")


    # save figures
    fig_rho.savefig("fig/" + filename + "_rho" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_enc.savefig("fig/" + filename + "_enc" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_sig.savefig("fig/" + filename + "_sig" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_Sig.savefig("fig/" + filename + "_Sig" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if kind > sphe:
        fig_ver.savefig("fig/" + filename + "_ver" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
        fig_sig2D.savefig("fig/" + filename + "_sigRz" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if outputPDF:
        fig_rho.savefig("fig/" + filename + "_rho" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_enc.savefig("fig/" + filename + "_enc" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_sig.savefig("fig/" + filename + "_sig" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_Sig.savefig("fig/" + filename + "_Sig" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        if kind > sphe:
            fig_ver.savefig("fig/" + filename + "_ver" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
            fig_sig2D.savefig("fig/" + filename + "_sigRz" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")

    if osipkov_merritt:
        fig_bet.savefig("fig/" + filename + "_bet" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
        if outputPDF:
            fig_bet.savefig("fig/" + filename + "_bet" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")

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
txtfile.close()


# read number of all component(s)
txtfile = open("doc/" + filename + ".summary.txt", "r")
unit = int(txtfile.readline())
line = txtfile.readline()
item = line.split("\t")
Nkind = int(item[0])
Nsphe = int(item[1])
txtfile.close()


# read analytic profile of all component(s) if possible
sphefile = "dat/" + filename + ".profile.h5"
add_ini_sphe = path.isfile(sphefile)
if add_ini_sphe:
    data_file = h5py.File(sphefile, "r")
    sphe_Nanal = data_file["/"].attrs["kinds"][0]
    sphe_Ndata = data_file["/data0/"].attrs["num"][0]
    sphe_rad = [0] * sphe_Nanal * sphe_Ndata
    sphe_rho = [0] * sphe_Nanal * sphe_Ndata
    sphe_enc = [0] * sphe_Nanal * sphe_Ndata
    sphe_Sig = [0] * sphe_Nanal * sphe_Ndata
    sphe_sig = [0] * sphe_Nanal * sphe_Ndata
    if osipkov_merritt:
        sphe_sgt = [0] * sphe_Nanal * sphe_Ndata
        sphe_bet = [0] * sphe_Nanal * sphe_Ndata
    else:
        sphe_sgt = [0]
        sphe_bet = [0]
    for kk in range(sphe_Nanal):
        folder = "data" + str(kk) + "/"
        sphe_rad[kk] = data_file[folder + "rad"]
        sphe_rho[kk] = data_file[folder + "rho"]
        sphe_enc[kk] = data_file[folder + "enc"]
        sphe_Sig[kk] = data_file[folder + "Sigma"]
        sphe_sig[kk] = data_file[folder + "sigma_r"]
        if osipkov_merritt:
            sphe_sgt[kk] = data_file[folder + "sigma_t"]
            sphe_bet[kk] = data_file[folder + "beta"]
    data_file.close()
else:
    sphe_rad = [0]
    sphe_rho = [0]
    sphe_enc = [0]
    sphe_Sig = [0]
    sphe_sig = [0]
    sphe_sgt = [0]
    sphe_bet = [0]


# diskfile = "dat/" + filename + ".disk.h5"
# add_ini_disk = path.isfile(diskfile)


cores = int(np.ceil(mp.cpu_count() / 2))
pool = mp.Pool(cores)
args = [(ii, Nkind, Nsphe, radmin, radmax, rhomin, rhomax, encmin, encmax, sigmin, sigmax, betmin, betmax, hormin, hormax, Sigmin, Sigmax, sigRmin, sigRmax, sigpmin, sigpmax, sigzmin, sigzmax, diskmin, diskmax, add_ini_sphe, sphe_rad, sphe_rho, sphe_enc, sphe_Sig, sphe_sig, sphe_sgt, sphe_bet) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
