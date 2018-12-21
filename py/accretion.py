import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import utils as utils


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", "--", ":", "-."]
Ncol, col = utils.set_color_palette_for_color_universal_design()


outputPDF = True


Nfile = 5
series = "m15"
# files = [series + "iso", series + "ra1_8", series + "ra1_4", series + "ra1_2", series + "ra1", series + "ra2"]
# lab = ["isotropic", r"$r_\mathrm{a} = r_\mathrm{h} / 8$", r"$r_\mathrm{a} = r_\mathrm{h} / 4$", r"$r_\mathrm{a} = r_\mathrm{h} / 2$", r"$r_\mathrm{a} = r_\mathrm{h}$", r"$r_\mathrm{a} = 2 r_\mathrm{h}$"]
files = [series + "iso", series + "ra1_4", series + "ra1_2", series + "ra1", series + "ra2"]
lab = ["isotropic", r"$r_\mathrm{a} = r_\mathrm{h} / 4$", r"$r_\mathrm{a} = r_\mathrm{h} / 2$", r"$r_\mathrm{a} = r_\mathrm{h}$", r"$r_\mathrm{a} = 2 r_\mathrm{h}$"]
component = ["DM halo", "bulge"]


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
plt.rcParams['font.size'] = 18
# plt.rcParams['font.size'] = 16

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


fig_menc = utils.set_figure(1, 1)
ax_menc = [0]
utils.locate_panels(fig_menc, ax_menc, 1, 1, True, True)
fig_frac = utils.set_figure(1, 1)
ax_frac = [0]
utils.locate_panels(fig_frac, ax_frac, 1, 1, True, True)
fig_mass = utils.set_figure(1, 1)
ax_mass = [0]
utils.locate_panels(fig_mass, ax_mass, 1, 1, True, True)

# mass_max = 0.0
for ii in range(Nfile):
    filename = "dat/" + files[ii] + ".losscone.h5"
    data = h5py.File(filename, "r")
    kind = data["/"].attrs["kind"][0]

    if ii == 0:
        nrad = data["/"].attrs["nrad"][0]
        rad = [0] * nrad
        enc = [0] * nrad
        rad = data["/rad"].value
        for kk in range(kind):
            enc = data["/Menc_" + str(kk)].value
            ax_menc[0].plot(rad, enc, linestyle = ls[kk % Nls], color = col[kk % Ncol], label = r"$M_\mathrm{orig}(r)$" + " (" + component[kk] + ")")

    for kk in range(kind):
        folder = "data" + str(kk) + "/"

        nrad = data[folder].attrs["nrad"][0]
        rad = [0] * nrad
        loss = [0] * nrad
        frac = [0] * nrad
        Msum = [0] * nrad

        rad = data[folder + "rad"].value
        loss = data[folder + "Mloc"].value
        frac = data[folder + "frac"].value
        Msum = data[folder + "Msum"].value

        jj = ii + kk * Nfile
        ax_menc[0].plot(rad, Msum, linestyle = ls[(jj + kind) % Nls], color = col[(jj + kind) % Ncol], label = r"$M_\mathrm{lost}(r)$" + " (" + lab[ii] + ", " + component[kk] + ")")

        ax_frac[0].plot(rad, frac, linestyle = ls[jj % Nls], color = col[jj % Ncol], label = lab[ii] + " (" + component[kk] + ")")

        ax_mass[0].plot(rad, loss, linestyle = ls[jj % Nls], color = col[jj % Ncol], label = lab[ii] + " (" + component[kk] + ")")
        # mass_max = max([mass_max, max(loss)])

    data.close()


ax_menc[0].set_xlabel(r"$r$ (kpc)")
ax_menc[0].set_ylabel(r"$M_\mathrm{enc}(r)$ ($M_\odot$)")
ax_frac[0].set_xlabel(r"$r$ (kpc)")
ax_frac[0].set_ylabel(r"Accreted fraction")
ax_mass[0].set_xlabel(r"$r$ (kpc)")
ax_mass[0].set_ylabel(r"$M(r)$ ($M_\odot$)")


# ax_mass[0].set_ylim(utils.scale_axis(mass_max * 3.1e-4, mass_max, True))


ax_menc[0].loglog()
ax_menc[0].grid()
ax_frac[0].loglog()
ax_frac[0].grid()
ax_mass[0].loglog()
ax_mass[0].grid()


# add legends
handles, labels = ax_menc[0].get_legend_handles_labels()
ax_menc[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')
handles, labels = ax_frac[0].get_legend_handles_labels()
ax_frac[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')
handles, labels = ax_mass[0].get_legend_handles_labels()
ax_mass[0].legend(handles, labels, numpoints = 1, handlelength = 2.5, loc = 'best')


# save figures
fig_menc.savefig("fig/" + series + "_accretion_menc" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
fig_frac.savefig("fig/" + series + "_accretion_frac" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
fig_mass.savefig("fig/" + series + "_accretion_mass" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
if outputPDF:
    fig_menc.savefig("fig/" + series + "_accretion_menc" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    fig_frac.savefig("fig/" + series + "_accretion_frac" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    fig_mass.savefig("fig/" + series + "_accretion_mass" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.close("all")
