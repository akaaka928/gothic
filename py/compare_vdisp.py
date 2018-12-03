import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import utils as utils


outputPDF = False


osipkov_merritt = False


Nfile = 2
files = ["hernquist", "hernquist_bh"]
tag = ["w/o BH", "w/ BH"]
Nkind = 1
lab = ["Hernquist"]
radmin, radmax = 1.0e-4, 15.0
sigmin, sigmax = 0.0, 4.5
osipkov_merritt = True


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", ":", "-.", "--"]
Ncol = 8
col = ["black", "red", "blue", "magenta", "green", "brown", "cyan", "yellow"]


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 14

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


fig = utils.set_figure(1, 1)
ax = [0]
utils.locate_panels(fig, ax, 1, 1, True, True)


for ii in range(Nfile):
    # read analytic profiles
    sphefile = "dat/" + files[ii] + ".profile.h5"
    data_file = h5py.File(sphefile, "r")
    sphe_Nanal = data_file["/"].attrs["kinds"][0]
    sphe_Ndata = data_file["/data0/"].attrs["num"][0]
    sphe_rad = [0] * sphe_Nanal * sphe_Ndata
    sphe_sig = [0] * sphe_Nanal * sphe_Ndata
    if osipkov_merritt:
        sphe_sgt = [0] * sphe_Nanal * sphe_Ndata
    else:
        sphe_sgt = [0]
    for kk in range(sphe_Nanal):
        folder = "data" + str(kk) + "/"
        sphe_rad[kk] = data_file[folder + "rad"].value
        sphe_sig[kk] = data_file[folder + "sigma_r"].value
        if osipkov_merritt:
            sphe_sgt[kk] = data_file[folder + "sigma_t"].value

    # read attributes
    length_unit = data_file["/"].attrs["length_astro_unit_name"][0]
    velocity_unit = data_file["/"].attrs["velocity_astro_unit_name"][0]
    data_file.close()


    for kk in range(Nkind):
        jj = kk + Nkind * ii
        if osipkov_merritt:
            ax[0].plot(sphe_rad[kk], sphe_sig[kk], linestyle = ls[(2 * jj    ) % Nls], color = col[jj % Ncol], label = r"$\sigma_r$" + " (" + lab[kk] + ", " + tag[ii] + ")")
            ax[0].plot(sphe_rad[kk], sphe_sgt[kk], linestyle = ls[(2 * jj + 1) % Nls], color = col[jj % Ncol], label = r"$\sigma_t$" + " (" + lab[kk] + ", " + tag[ii] + ")")
        else:
            ax[0].plot(sphe_rad[kk], sphe_sig[kk], linestyle = ls[jj % Nls], color = col[jj % Ncol])


ax[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode("UTF-8")))
if osipkov_merritt:
    ax[0].set_ylabel(r"$\sigma$ ({:<})".format(velocity_unit.decode("UTF-8")))
else:
    ax[0].set_ylabel(r"$\sigma_r$ ({:<})".format(velocity_unit.decode("UTF-8")))

ax[0].set_xlim(utils.scale_axis(radmin, radmax,  True))
ax[0].set_ylim(utils.scale_axis(sigmin, sigmax, False))


ax[0].semilogx()
ax[0].grid()


# add legends
handles, labels = ax[0].get_legend_handles_labels()
ax[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')


# save figures
fig.savefig("fig/" + files[0] + "_compsig" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
if outputPDF:
    fig.savefig("fig/" + files[0] + "_compsig" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.close("all")
