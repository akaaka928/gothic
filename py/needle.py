import h5py
import numpy as np

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import utils as utils


nxpanel, nypanel = 2, 1
filename = "cb17"
col = ["black", "red", "blue", "magenta"]
khead = 1


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 24


fig = utils.set_figure(nxpanel, nypanel)
ax = [0] * nxpanel * nypanel
utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)


# set plot range
# xmin, xmax = -25.0, 25.0
# zmin, zmax = -25.0, 25.0
xmin, xmax = -15.0, 15.0
zmin, zmax = -15.0, 15.0
# xmin, xmax = -10.0, 10.0
# zmin, zmax = -10.0, 10.0

for ii in range(nxpanel):
    idx = ii
    jj = 0
    # pick up an appropriate snapshot
    if (ii & 1) == 0:
        model = "1kpc"
    else:
        model = "needle/1kpc"
    input_file = model + "/dat/" + filename + ".split" + "000" + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    length_unit = h5file["/"].attrs["length_astro_unit_name"]
    time_unit = h5file["/"].attrs["time_astro_unit_name"]
    mass_unit = h5file["/"].attrs["mass_astro_unit_name"]
    time = h5file["/"].attrs["time"]
    kind = int(h5file["/"].attrs["kinds"])

    for kk in range(khead, kind):
        # read particle position and mass
        folder = "data" + str(kk) + "/"
        position = h5file[folder + "position"].value

        # data preparation
        px = position[:, 0]
        py = position[:, 1]
        pz = position[:, 2]

        # plot the data
        # ax[idx].plot(px, pz, ",", color = col[kk - khead], rasterized = True)
        ax[idx].plot(px, pz, ".", markersize = 1, color = col[kk - khead], rasterized = True)


    # close the HDF5 file
    h5file.close()

    # set plot range
    ax[idx].set_xlim([xmin, xmax])
    ax[idx].set_ylim([zmin, zmax])

    ax[idx].set_xticks([-10, -5, 0, 5, 10])
    ax[idx].set_yticks([-10, -5, 0, 5, 10])
    # ax[idx].set_xticks([-9, -6, -3, 0, 3, 6, 9])
    # ax[idx].set_yticks([-9, -6, -3, 0, 3, 6, 9])
    # ax[idx].grid()
    ax[idx].tick_params(axis = "both", direction = "in", color = "black", bottom = "on", top = "on", left = "on", right = "on")

    # set label
    if jj == 0:
        ax[idx].set_xlabel(r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
    if ii == 0:
        ax[idx].set_ylabel(r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))

    # set caption
    caption  = "(" + "{:^c}".format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ")"
    ax[idx].text(xmin + 1, zmax - 2, caption, fontsize=20)


utils.set_shared_xlabel(ax, r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
utils.set_shared_ylabel(ax, r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))


# plt.show()
plt.savefig("needle.png", format = "png", dpi = 300, bbox_inches = "tight")
plt.savefig("needle.pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
