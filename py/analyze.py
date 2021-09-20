import numpy as np
import math

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import utils as utils


Nmass = 7
num  = [6, 6, 6, 6, 5, 5, 5]
head = [0, 6, 12, 18, 24, 29, 34]
lab  = [r"$M_{200} = 10^9 M_\odot$", r"$M_{200} = 10^{10} M_\odot$", r"$M_{200} = 10^{11} M_\odot$", r"$M_{200} = 10^{12} M_\odot$", r"$M_{200} = 10^{13} M_\odot$", r"$M_{200} = 10^{14} M_\odot$", r"$M_{200} = 10^{15} M_\odot$"]


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", "--", ":", "-."]
Ncol, col = utils.set_color_palette_for_color_universal_design()


outputPDF = True


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath,siunitx}"


# set font size
plt.rcParams['font.size'] = 28

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


nxpanel, nypanel = 1, 2
fig_mass_1ra = utils.set_figure(nxpanel, nypanel)
ax_mass_1ra = [0] * nxpanel * nypanel
utils.locate_panels(fig_mass_1ra, ax_mass_1ra, nxpanel, nypanel, True, True)
fig_mass_2ra = utils.set_figure(nxpanel, nypanel)
ax_mass_2ra = [0] * nxpanel * nypanel
utils.locate_panels(fig_mass_2ra, ax_mass_2ra, nxpanel, nypanel, True, True)
fig_size_1ra = utils.set_figure(nxpanel, nypanel)
ax_size_1ra = [0] * nxpanel * nypanel
utils.locate_panels(fig_size_1ra, ax_size_1ra, nxpanel, nypanel, True, True)
fig_size_2ra = utils.set_figure(nxpanel, nypanel)
ax_size_2ra = [0] * nxpanel * nypanel
utils.locate_panels(fig_size_2ra, ax_size_2ra, nxpanel, nypanel, True, True)


# read dataset
data = np.genfromtxt("./summary.txt", comments = "#", delimiter = "\t")
# print(data)
M200  = data[:,  0]
M_bh  = data[:,  1]
rs    = data[:,  2]
rhalf = data[:,  3]
r200  = data[:,  4]
ra    = data[:,  5]
Mget  = data[:,  6]
Mhalf = data[:,  7]
rmid  = data[:,  8]
tmid  = data[:,  9]
Mdot  = data[:, 10]
f1ra  = data[:, 11]
f2ra  = data[:, 12]
f3ra  = data[:, 13]

Mget_spec = Mget / M_bh
Mdot_spec = Mdot / M_bh

for ii in range(Nmass):
    jj = Nmass - ii - 1
    ax_mass_1ra[1].plot(f1ra[head[ii]:head[ii]+num[ii]], Mget_spec[head[ii]:head[ii]+num[ii]], linestyle = ls[jj % Nls], marker = pt[jj % Npt], color = col[jj % Ncol], label = lab[ii])
    ax_mass_1ra[0].plot(f1ra[head[ii]:head[ii]+num[ii]], Mdot_spec[head[ii]:head[ii]+num[ii]], linestyle = ls[jj % Nls], marker = pt[jj % Npt], color = col[jj % Ncol], label = lab[ii])
    ax_mass_2ra[1].plot(f2ra[head[ii]:head[ii]+num[ii]], Mget_spec[head[ii]:head[ii]+num[ii]], linestyle = ls[jj % Nls], marker = pt[jj % Npt], color = col[jj % Ncol], label = lab[ii])
    ax_mass_2ra[0].plot(f2ra[head[ii]:head[ii]+num[ii]], Mdot_spec[head[ii]:head[ii]+num[ii]], linestyle = ls[jj % Nls], marker = pt[jj % Npt], color = col[jj % Ncol], label = lab[ii])
    ax_size_1ra[1].plot(f1ra[head[ii]:head[ii]+num[ii]],      rmid[head[ii]:head[ii]+num[ii]], linestyle = ls[jj % Nls], marker = pt[jj % Npt], color = col[jj % Ncol], label = lab[ii])
    ax_size_1ra[0].plot(f1ra[head[ii]:head[ii]+num[ii]],      tmid[head[ii]:head[ii]+num[ii]], linestyle = ls[jj % Nls], marker = pt[jj % Npt], color = col[jj % Ncol], label = lab[ii])
    ax_size_2ra[1].plot(f2ra[head[ii]:head[ii]+num[ii]],      rmid[head[ii]:head[ii]+num[ii]], linestyle = ls[jj % Nls], marker = pt[jj % Npt], color = col[jj % Ncol], label = lab[ii])
    ax_size_2ra[0].plot(f2ra[head[ii]:head[ii]+num[ii]],      tmid[head[ii]:head[ii]+num[ii]], linestyle = ls[jj % Nls], marker = pt[jj % Npt], color = col[jj % Ncol], label = lab[ii])


for ii in range(nxpanel):
    ax_mass_1ra[ii * nypanel].set_xlabel(r"$M(r \leq r_\mathrm{a}) / M_{200}$")
    ax_mass_2ra[ii * nypanel].set_xlabel(r"$M(r \leq 2 r_\mathrm{a}) / M_{200}$")
    ax_size_1ra[ii * nypanel].set_xlabel(r"$M(r \leq r_\mathrm{a}) / M_{200}$")
    ax_size_2ra[ii * nypanel].set_xlabel(r"$M(r \leq 2 r_\mathrm{a}) / M_{200}$")
    for jj in range(nypanel):
        idx = jj + ii * nypanel
        ax_mass_1ra[idx].grid()
        ax_mass_1ra[idx].set_xscale("log")
        ax_mass_1ra[idx].set_yscale("log")
        ax_mass_2ra[idx].grid()
        ax_mass_2ra[idx].set_xscale("log")
        ax_mass_2ra[idx].set_yscale("log")
        ax_size_1ra[idx].grid()
        ax_size_1ra[idx].set_xscale("log")
        ax_size_1ra[idx].set_yscale("log")
        ax_size_2ra[idx].grid()
        ax_size_2ra[idx].set_xscale("log")
        ax_size_2ra[idx].set_yscale("log")

ax_mass_1ra[1].set_ylabel(r"$M_\mathrm{accreted} / M_\mathrm{BH}$")
ax_mass_1ra[0].set_ylabel(r"$\dot{M} / M_\mathrm{BH} \, (\si{yr^{-1}})$")
ax_mass_2ra[1].set_ylabel(r"$M_\mathrm{accreted} / M_\mathrm{BH}$")
ax_mass_2ra[0].set_ylabel(r"$\dot{M} / M_\mathrm{BH} \, (\si{yr^{-1}})$")

ax_size_1ra[1].set_ylabel(r"$r_{1/2}~(\si{kpc})$")
ax_size_1ra[0].set_ylabel(r"$t_{1/2}~(\si{Myr})$")
ax_size_2ra[1].set_ylabel(r"$r_{1/2}~(\si{kpc})$")
ax_size_2ra[0].set_ylabel(r"$t_{1/2}~(\si{Myr})$")

# add legends
handles, labels = ax_mass_1ra[1].get_legend_handles_labels()
ax_mass_1ra[1].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best', fontsize = 18)
handles, labels = ax_mass_2ra[1].get_legend_handles_labels()
ax_mass_2ra[1].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best', fontsize = 18)
handles, labels = ax_size_1ra[1].get_legend_handles_labels()
ax_size_1ra[1].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best', fontsize = 18)
handles, labels = ax_size_2ra[1].get_legend_handles_labels()
ax_size_2ra[1].legend(handles[::-1], labels[::-1], numpoints = 1, handlelength = 2.5, loc = 'best', fontsize = 18)


# save figures
fig_mass_1ra.savefig("fig/losscone_mass_1ra" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
fig_mass_2ra.savefig("fig/losscone_mass_2ra" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
fig_size_1ra.savefig("fig/losscone_size_1ra" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
fig_size_2ra.savefig("fig/losscone_size_2ra" + ".png", format = "png", dpi =  96, bbox_inches = "tight")
if outputPDF:
    fig_mass_1ra.savefig("fig/losscone_mass_1ra" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    fig_mass_2ra.savefig("fig/losscone_mass_2ra" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    fig_size_1ra.savefig("fig/losscone_size_1ra" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    fig_size_2ra.savefig("fig/losscone_size_2ra" + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")


plt.close("all")
