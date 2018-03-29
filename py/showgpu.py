import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import os.path as path
import multiprocessing as mp

import utils as utils


outputPDF = True
multiGPUs = True

# filename = "gss"
# init = 0
# last = 6

# filename = "k17disk"
# init = 0
# last = 399
# multiGPUs = False

filename = "halocusp_run"
init = 0
last = 140
multiGPUs = False


pt = ["o", "s", "^", "D"]
ls = ["-", ":", "-.", "--"]
col = ["black", "red", "blue", "magenta"]


def draw_figure(fileid):
    # read snapshot
    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".snp" + snapshot + ".h5"
    h5file = h5py.File(input_file, "r")
    steps = h5file["snp"].attrs["steps"][0]
    record = h5file["GPUinfo/record"].value
    # print(record.shape)
    # print(record.dtype)
    h5file.close()

    if multiGPUs:
        Ngpus = record.shape[0]
        Nstep = record.shape[1]
    else:
        Ngpus = 1
        Nstep = record.shape[0]

    elapsed     = np.zeros((Ngpus, Nstep))
    devClock    = np.zeros((Ngpus, Nstep))
    temperature = np.zeros((Ngpus, Nstep))
    power       = np.zeros((Ngpus, Nstep))
    grpNum      = np.zeros((Ngpus, Nstep))

    if multiGPUs:
        for ii in range(Ngpus):
            for jj in range(Nstep):
                elapsed[ii][jj]     = record[ii][jj][0]
                devClock[ii][jj]    = record[ii][jj][1]
                temperature[ii][jj] = record[ii][jj][2]
                power[ii][jj]       = record[ii][jj][3]
                grpNum[ii][jj]      = record[ii][jj][4]
    else:
        for jj in range(Nstep):
            elapsed[0][jj]     = record[jj][0]
            devClock[0][jj]    = record[jj][1]
            temperature[0][jj] = record[jj][2]
            power[0][jj]       = record[jj][3]
            grpNum[0][jj]      = record[jj][4]


    step = [0] * Nstep
    for ii in range(Nstep):
        step[ii] = steps - (Nstep - 1 - ii)

    # plot the data
    nxpanel = 1
    nypanel = 1
    fig_time = utils.set_figure(nxpanel, nypanel)
    fig_clck = utils.set_figure(nxpanel, nypanel)
    fig_temp = utils.set_figure(nxpanel, nypanel)
    fig_powr = utils.set_figure(nxpanel, nypanel)
    fig_gnum = utils.set_figure(nxpanel, nypanel)
    ax_time = [0] * nxpanel * nypanel
    ax_clck = [0] * nxpanel * nypanel
    ax_temp = [0] * nxpanel * nypanel
    ax_powr = [0] * nxpanel * nypanel
    ax_gnum = [0] * nxpanel * nypanel
    utils.locate_panels(fig_time, ax_time, nxpanel, nypanel, True, True)
    utils.locate_panels(fig_clck, ax_clck, nxpanel, nypanel, True, True)
    utils.locate_panels(fig_temp, ax_temp, nxpanel, nypanel, True, True)
    utils.locate_panels(fig_powr, ax_powr, nxpanel, nypanel, True, True)
    utils.locate_panels(fig_gnum, ax_gnum, nxpanel, nypanel, True, True)

    for ii in range(Ngpus):
        ax_time[0].plot(step,     elapsed[ii], "x", linestyle = "-", linewidth = 0.5)
        ax_clck[0].plot(step,    devClock[ii], "x", linestyle = "-", linewidth = 0.5)
        ax_temp[0].plot(step, temperature[ii], "x", linestyle = "-", linewidth = 0.5)
        ax_powr[0].plot(step,       power[ii], "x", linestyle = "-", linewidth = 0.5)
        ax_gnum[0].plot(step,      grpNum[ii], "x", linestyle = "-", linewidth = 0.5)

    ax_time[0].set_xlabel(r"time step")
    ax_clck[0].set_xlabel(r"time step")
    ax_temp[0].set_xlabel(r"time step")
    ax_powr[0].set_xlabel(r"time step")
    ax_gnum[0].set_xlabel(r"time step")

    ax_time[0].set_ylabel(r"execution time (\si{\second})")
    ax_clck[0].set_ylabel(r"clock frequency (\si{\mega\hertz})")
    ax_temp[0].set_ylabel(r"device temperature (\si{\degreeCelsius})")
    ax_powr[0].set_ylabel(r"power (\si{\milli\watt})")
    ax_gnum[0].set_ylabel(r"number of particle groups")

    ax_time[0].semilogy()
    ax_gnum[0].semilogy()

    ax_time[0].grid()
    ax_clck[0].grid()
    ax_temp[0].grid()
    ax_powr[0].grid()
    ax_gnum[0].grid()


    # save figures
    figname = "fig/" + filename + "_map" + snapshot
    fig_time.savefig("fig/" + filename + "_time" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_clck.savefig("fig/" + filename + "_clck" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_temp.savefig("fig/" + filename + "_temp" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_powr.savefig("fig/" + filename + "_powr" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    fig_gnum.savefig("fig/" + filename + "_gnum" + snapshot + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if outputPDF:
        fig_time.savefig("fig/" + filename + "_time" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_clck.savefig("fig/" + filename + "_clck" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_temp.savefig("fig/" + filename + "_temp" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_powr.savefig("fig/" + filename + "_powr" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
        fig_gnum.savefig("fig/" + filename + "_gnum" + snapshot + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.close("all")



# def wrapper(argv):
#     return draw_figure(*argv)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# set font size
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 14

plt.rcParams['text.latex.preamble'] = [r"\usepackage{siunitx}"]

# cores = mp.cpu_count()
cores = int(np.ceil(mp.cpu_count() / 2))
pool = mp.Pool(cores)
# args = [(ii, Nkind) for ii in range(init, last + 1, 1)]
# pool.map(wrapper, args)
pool.map(draw_figure, range(init, last + 1, 1))
pool.close()
