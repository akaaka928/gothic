import sys, os
import array
import numpy as np
import matplotlib
# matplotlib.use("tkagg")
matplotlib.use("agg")
import matplotlib.pyplot as plt

import utils as utils

# filename = "cb17"
# kind = 2
filename = "m31"
kind = 1

xmax, zmax = 21.0, 1.2


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
# plt.rcParams['font.size'] = 20
plt.rcParams['font.size'] = 28

# set number of panels
nxpanel, nypanel = 1, 1
fig = utils.set_figure(nxpanel, nypanel)
ax = [0] * nxpanel * nypanel
utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)


lines = [] * 2
for ii in range(kind):
    for jj in range(2):
        # read particle data
        input_file = "dat/" + filename + ".disk" + str(ii) + ".snp0" + str(jj * 4) + "0.dat"

        fsize = os.path.getsize(input_file)
        num = int(np.rint(fsize / (4 * 2))) # 2 data (radius and height) in FP32 (4 bytes)
        rad = [0] * num
        rms = [0] * num

        fin = open(input_file, "rb")
        arr = array.array("f")
        arr.fromfile(fin, num * 2)
        fin.close()

        for kk in range(num):
            rad[kk] = arr[      kk]
            rms[kk] = arr[num + kk]


        # plot the data
        if ii == 0:
            col = "black"
        else:
            col = "red"

        if jj == 0:
            sty = ":"
            wid = 1
        else:
            sty = "-"
            wid = 2

        if ii == 0:
            lines += ax[0].plot(rad, rms, color = col, linestyle = sty, linewidth = wid)
        else:
            ax[0].plot(rad, rms, color = col, linestyle = sty, linewidth = wid)

# set plot range
ax[0].set_xlim([0, xmax])
ax[0].set_ylim([0, zmax])

ax[0].grid()
ax[0].tick_params(axis = "both", direction = "in", color = "black", bottom = "on", top = "on", left = "on", right = "on")


ax[0].text(10.5, 0.42, "thin disc", color = "black")
ax[0].text(10.5, 0.90, "thick disc", color = "black")


# set label
ax[0].set_xlabel(r"$R$ (kpc)")
ax[0].set_ylabel(r"$\mathrm{RMS}(z)$ (kpc)")

labels = [r"$t = 1$~Gyr", r"$t = 0$~Gyr"]
ax[0].legend(lines[::-1], labels)

figname = "fig/" + filename + "_rms"
plt.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
plt.savefig(figname + ".png", format = "png", dpi = 300, bbox_inches = "tight")
