import h5py
import numpy as np

import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt


def locate_panels(ax, nx, ny, share_xaxis, share_yaxis):
    margin = 0.12
    if (share_xaxis == False) or (share_yaxis == False):
        margin = 0.15

    xmin, xmax = margin, 1.0 - margin
    ymin, ymax = margin, 1.0 - margin
    xbin = (xmax - xmin) / nx
    ybin = (ymax - ymin) / ny
    xmargin, ymargin = 0, 0

    if share_yaxis == False:
        xmin = 0.0
        xbin = 1.0 / nx
        xmargin = xbin * margin

    if share_xaxis == False:
        ymin = 0.0
        ybin = 1.0 / ny
        ymargin = ybin * margin

    for ii in range(nx):
        xl = xmin + ii * xbin + xmargin

        for jj in range(ny):
            yl = ymin + jj * ybin + ymargin
            kk = ii * ny + jj
            ax[kk] = fig.add_axes((xl, yl, xbin - 2 * xmargin, ybin - 2 * ymargin))

            if share_xaxis == True:
                ax[kk].tick_params(labelbottom = "off")
                if jj == 0:
                    ax[kk].tick_params(labelbottom = "on")

            if share_yaxis == True:
                ax[kk].tick_params(labelleft = "off")
                if ii == 0:
                    ax[kk].tick_params(labelleft = "on")


def set_shared_xlabel(ax, xlabel):
    fig = ax[-1].figure
    fig.canvas.draw()

    # get the corner for all plots
    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax:
        at.set_xlabel("") # remove existing xlabels

        xl, xr = at.get_position().x0, at.get_position().x1

        bboxes, _ = at.xaxis.get_ticklabel_extents(fig.canvas.renderer)
        bboxes = bboxes.inverse_transformed(fig.transFigure)
        yb, yt = bboxes.y0, bboxes.y1

        if x0 > xl:
            x0 = xl
        if x1 < xr:
            x1 = xr
        if y0 > yb:
            y0 = yb
        if y1 < yt:
            y1 = yt

    # set position of label
    ax[-1].set_xlabel(xlabel)
    ax[-1].xaxis.set_label_coords((x0 + x1) / 2, (y0 + y1) / 2, transform = fig.transFigure)


def set_shared_ylabel(ax, ylabel):
    fig = ax[-1].figure
    fig.canvas.draw()

    # get the corner for all plots
    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax:
        at.set_ylabel("") # remove existing ylabels

        yb, yt = at.get_position().y0, at.get_position().y1

        bboxes, _ = at.yaxis.get_ticklabel_extents(fig.canvas.renderer)
        bboxes = bboxes.inverse_transformed(fig.transFigure)
        xl, xr = bboxes.x0, bboxes.x1

        if x0 > xl:
            x0 = xl
        if x1 < xr:
            x1 = xr
        if y0 > yb:
            y0 = yb
        if y1 < yt:
            y1 = yt

    # set position of label
    ax[-1].set_ylabel(ylabel)
    ax[-1].yaxis.set_label_coords((x0 + x1) / 2, (y0 + y1) / 2, transform = fig.transFigure)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 14


# set number of panels
nxpanel, nypanel = 4, 4
ax = [0] * nxpanel * nypanel

# set figure size and its aspect ratio
Lx = 10
if nxpanel > 2 * nypanel:
    Lx *= 2
Ly = (Lx / nxpanel) * nypanel
fig = plt.figure(figsize = (Lx, Ly))

# set location of panels
locate_panels(ax, nxpanel, nypanel, True, True)

# set plot range
xmin, xmax = -25.0, 25.0
zmin, zmax = -25.0, 25.0

for ii in range(nxpanel):
    # pick up an appropriate snapshot
    if (ii & 1) == 0:
        snapshot = "000"
    else:
        snapshot = "040"
    for jj in range(nypanel):
        idx = ii * nypanel + jj
        kk = jj * (nxpanel >> 1) + ((nxpanel - 1 - ii) >> 1)
        param = ((nxpanel >> 1) * nypanel - 1 - kk) + 1
        model = str(param) + "kpc"
        input_file = model + "/dat/ltg.split" + snapshot + ".h5"

        # read snapshot
        h5file = h5py.File(input_file, "r")

        # read attributes
        length_unit = h5file["/"].attrs["length_astro_unit_name"]
        time_unit = h5file["/"].attrs["time_astro_unit_name"]
        mass_unit = h5file["/"].attrs["mass_astro_unit_name"]
        time = h5file["/"].attrs["time"]

        # # confirmation of attributes
        # print("{:.1f} {:<}".format(time[0], time_unit[0].decode('UTF-8')))
        # print("unit of the length is {:<}".format(length_unit[0].decode('UTF-8')))
        # print("unit of the mass is {:<}".format(mass_unit[0].decode('UTF-8')))

        # read particle position and mass
        folder = "data2/"
        position = h5file[folder + "position"].value
        mass = h5file[folder + "mass"].value

        # close the HDF5 file
        h5file.close()

        # data preparation
        px = position[:, 0]
        py = position[:, 1]
        pz = position[:, 2]

        # plot the data
        ax[idx].plot(px, pz, ",", color = "black", rasterized = True)

        # set plot range
        ax[idx].set_xlim([xmin, xmax])
        ax[idx].set_ylim([zmin, zmax])

        ax[idx].set_xticks([-20, -10, 0, 10, 20])
        ax[idx].set_yticks([-20, -10, 0, 10, 20])
        # ax[idx].grid()
        ax[idx].tick_params(axis = "both", direction = "in", color = "black", bottom = "on", top = "on", left = "on", right = "on")

        # set label
        if jj == 0:
            ax[idx].set_xlabel(r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
        if ii == 0:
            ax[idx].set_ylabel(r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))

        # set caption
        caption  = "(" + "{:^c}".format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ")"
        caption += " " + r"$z_\mathrm{d}=$" + r"${:.0f}$ {:<}".format(param, length_unit[0].decode('UTF-8'))
        caption += ", $t =$" + r"${:.0f}$ {:<}".format(time[0] / 1000, "Gyr")
        ax[idx].text(xmin + 2, zmax - 7, caption, fontsize=12)


set_shared_xlabel(ax, r"$x$ ({:<})".format(length_unit[0].decode('UTF-8')))
set_shared_ylabel(ax, r"$z$ ({:<})".format(length_unit[0].decode('UTF-8')))


# plt.show()
plt.savefig("dot.png", format = "png", dpi = 300, bbox_inches = "tight")
plt.savefig("dot.pdf", format = "pdf", dpi = 300, bbox_inches = "tight") # too slow
# plt.savefig("dot.svg", format = "svg", dpi = 300, bbox_inches = "tight") # too slow
