import numpy as np
import h5py

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt

import utils as utils


# specify plot target
# filename = "cb17"
filename = "m31"

# conversion factors of unit system (taken from doc/unit.txt)
newton_com = 4.498466e+00
mass2astro = 1.000000e+08
length2astro = 1.000000e+00
velocity2astro = 9.777922e+00


# set plot range
# rmin, rmax = 1.0e-1, 2.0e+2
# vmin, vmax = 0.0, 200.0
rmin, rmax = 1.0e-1, 4.0e+2
vmin, vmax = 0.0, 300.0


# set number of panels
nxpanel, nypanel = 1, 1


col = ["black", "blue", "magenta", "green"]
ls  = ["-", "-.", ":", "--"]
# lab = ["dark matter halo", "bulge", "thick disc", "thin disc"]
lab = ["dark matter halo", "bulge", "disk"]


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 24
# plt.rcParams['font.size'] = 28


# read analytic profile of all component(s)
target = "dat/" + filename + ".profile.h5"
data_file = h5py.File(target, "r")
length_unit = data_file["/"].attrs["length_astro_unit_name"][0]
velocity_unit = data_file["/"].attrs["velocity_astro_unit_name"][0]
Nkind = data_file["/"].attrs["kinds"][0]
Ndata = data_file["/data0/"].attrs["num"][0]
rad = np.empty(Nkind * Ndata)
enc = np.empty(Nkind * Ndata)
for kk in range(Nkind):
    folder = "data" + str(kk) + "/"
    rr = data_file[folder + "rad"].value
    MM = data_file[folder + "enc"].value
    for ii in range(Ndata):
        rad[kk * Ndata + ii] = rr[ii]
        enc[kk * Ndata + ii] = MM[ii]
data_file.close()

conversion = np.sqrt(newton_com * length2astro / mass2astro) * velocity2astro

vel = np.empty(Nkind * Ndata)
vel = np.sqrt(enc / rad) * conversion


radius = [0] * Ndata
vcirc = [0] * Ndata
for ii in range(Ndata):
    radius[ii] = rad[ii]
    mm = 0
    for kk in range(Nkind):
        mm += enc[kk * Ndata + ii]
    vcirc[ii] = np.sqrt(mm / radius[ii]) * conversion


fig = utils.set_figure(nxpanel, nypanel)
ax = [0] * nxpanel * nypanel
utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

for ll in range(Nkind):
    kk = Nkind - 1 - ll
    ax[0].plot(rad[kk * Ndata : (kk + 1) * Ndata - 1], vel[kk * Ndata : (kk + 1) * Ndata - 1], linestyle = ls[kk], color = col[kk], label = lab[kk])
ax[0].plot(radius, vcirc, linestyle = "-", color = "red", label = r"total")

# set plot range
ax[0].set_xlim([rmin, rmax])
ax[0].set_ylim([vmin, vmax])
ax[0].semilogx()
ax[0].grid()

# set label
ax[0].set_xlabel(r"$r$ ({:<})".format(length_unit.decode('UTF-8')))
# ax[0].set_ylabel(r"$v_\mathrm{rot}$" + " ({:<})".format(velocity_unit.decode('UTF-8')))
ax[0].set_ylabel(r"$v_\mathrm{rot}$ (km s$^{-1}$)")

# set legend
handles, labels = ax[0].get_legend_handles_labels()
ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, prop = {'size' : 20}, loc = 'best')
# ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best')

# output figure
figname = "fig/" + filename + "_vrot"
plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")
