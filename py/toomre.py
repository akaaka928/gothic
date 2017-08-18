import numpy as np
import h5py

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt

import utils as utils


# specify plot target
filename = "cb17"

# set plot range
Rmin, Rmax = 1.0e-1, 50.0
Qmin, Qmax = 1.0e-1, 1.0e+2


# set number of panels
nxpanel, nypanel = 1, 1


# col = ["black", "red", "blue", "magenta", "green"]
# ls  = ["-", "-.", ":", "--"]
col = ["black", "black"]
ls  = ["-", ":"]
lab = ["thick disc", "thin disc"]


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
plt.rcParams['font.size'] = 28


# read analytic profile of all component(s)
target = "dat/" + filename + ".disk.h5"
data_file = h5py.File(target, "r")
length_unit = data_file["/"].attrs["length_astro_unit_name"][0]
Nkind = data_file["/"].attrs["kinds"][0]
Ndata = data_file["/data0/1D_data/"].attrs["num"][0]
RR = np.empty(Nkind * Ndata)
QQ = np.empty(Nkind * Ndata)
for kk in range(Nkind):
    folder = "data" + str(kk) + "/1D_data/"
    rad = data_file[folder + "radius"].value
    val = data_file[folder + "Toomre's Q"].value
    for ii in range(Ndata):
        RR[kk * Ndata + ii] = rad[ii]
        QQ[kk * Ndata + ii] = val[ii]
data_file.close()



fig = utils.set_figure(nxpanel, nypanel)
ax = [0] * nxpanel * nypanel
utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

for ll in range(Nkind):
    kk = Nkind - 1 - ll
    ax[0].plot(RR[kk * Ndata : (kk + 1) * Ndata - 1], QQ[kk * Ndata : (kk + 1) * Ndata - 1], linestyle = ls[kk], color = col[kk], label = lab[kk])

# set plot range
ax[0].set_xlim([Rmin, Rmax])
ax[0].set_ylim([Qmin, Qmax])
ax[0].loglog()
ax[0].grid()

# set label
ax[0].set_xlabel(r"$R$ ({:<})".format(length_unit.decode('UTF-8')))
# ax[0].set_ylabel(r"$v_\mathrm{rot}$" + " ({:<})".format(velocity_unit.decode('UTF-8')))
ax[0].set_ylabel(r"$Q$")

# set legend
handles, labels = ax[0].get_legend_handles_labels()
# ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, prop = {'size' : 14}, loc = 'best')
# ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, prop = {'size' : 14}, loc = 'upper left')
ax[0].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'upper left')

# output figure
figname = "fig/" + filename + "_Q"
plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")
