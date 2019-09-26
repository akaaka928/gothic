import numpy as np
import math
import h5py

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.cm import get_cmap # for setting facecolor in imshow
from matplotlib.colors import LogNorm # for logarithmic plot in imshow

import utils as utils


from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


from argparse import ArgumentParser
def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-f', '--files', type=str,
                           help='Name of the target files')
    argparser.add_argument('-c', '--continued',
                           action='store_true',
                           help='Visualize continued simulations')
    argparser.add_argument('-p', '--pdf',
                           action='store_true',
                           help='Whether to output figures in PDF format')
    return argparser.parse_args()

# catch optional input arguments
args = get_option()
filename = args.files
continued = args.continued
outputPDF = args.pdf


fs_base = 16
tl_base = 6.0
tw_base = 1.0


pt = ['o', 's', '^', 'D']
ls = ['-', ':', '-.', '--']
Ncol, col = utils.set_color_palette_for_color_universal_design()


time_offset = 0.0
if continued:
    time_offset = -775.0
    # time_offset = -4875.0


Nskip = 2
skip = [1, 3]
Nmbh = 0
init = 0
last = 264
# last = 0
lab = ['Leading group', 'Trailing group']
fvmin, fvmax = 1.0e+2, 1.0e+5
fEmin, fEmax = 1.0e-3, 1.0e+0



# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# use packages
plt.rcParams['text.latex.preamble'] = [r'\usepackage{physics}']
plt.rcParams['text.latex.preamble'] = [r'\usepackage{siunitx}']

# set font size
plt.rcParams['font.size'] = fs_base

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# read number of all component(s)
if rank == 0:
    with open('doc/' + filename + '.summary.txt', 'r') as txtfile:
        unit = int(txtfile.readline())
        line = txtfile.readline()
        item = line.split('\t')
        kind = int(item[0])
        # sphe = int(item[1])
else:
    kind = None
kind = comm.bcast(kind, root = 0)


# cmap_Sigma = 'plasma'
# cmap_depth = 'cividis'
# cmap_depth = 'viridis'
# cmap_phase = 'inferno'
cmap_vel = 'inferno'
cmap_ene = 'plasma'


for fileid in range(init + rank, last + 1, size):
    snapshot = '{:03d}'.format(fileid)
    input_file = 'dat/' + filename + '.m31ene' + snapshot + '.h5'

    # read snapshot
    with h5py.File(input_file, 'r') as h5file:
        # read attributes
        time = h5file['/'].attrs['time'][0] + time_offset

        # kind = h5file['/'].attrs['kinds'][0]
        nvr = h5file['/'].attrs['nvr'][0]
        nvt = h5file['/'].attrs['nvt'][0]
        nE = h5file['/'].attrs['nE'][0]
        nJ = h5file['/'].attrs['nJ'][0]

        # memory allocation for (vr, vt)-maps
        vrvt = np.zeros((kind + 1, nvr, nvt))
        vr = np.zeros((kind + 1, nvr + 1))
        vt = np.zeros((kind + 1, nvt + 1))

        # memory allocation for (E, J)-maps
        EJ = np.zeros((kind + 1, nJ, nE))
        EE = np.zeros((kind + 1, nE + 1))
        JJ = np.zeros((kind + 1, nJ + 1))

        nxpanel = 1
        nypanel = 0
        Ndata = 0
        for ii in range(kind):
            if (Nskip == 0) or (ii not in skip):
                # read maps and axes
                folder = 'map' + str(ii) + '/'
                vrvt[nypanel] = h5file[folder + 'vrvt']
                vr[nypanel] = h5file[folder + 'vr']
                vt[nypanel] = h5file[folder + 'vt']
                EJ[nypanel] = h5file[folder + 'EJ']
                EE[nypanel] = h5file[folder + 'E']
                JJ[nypanel] = h5file[folder + 'J']

                nypanel += 1
                Ndata += 1
        # close the HDF5 file


    vrvt = np.maximum(vrvt, fvmin)
    vrvt = np.minimum(vrvt, fvmax)
    EJ = np.minimum(EJ, fEmax)
    EJ = np.maximum(EJ, fEmin)


    # plot the data
    fig_v = utils.set_figure(nxpanel, nypanel)
    ax_v = [0] * nxpanel * nypanel
    utils.locate_panels(fig_v, ax_v, nxpanel, nypanel, True, True)
    fig_E = utils.set_figure(nxpanel, nypanel)
    ax_E = [0] * nxpanel * nypanel
    utils.locate_panels(fig_E, ax_E, nxpanel, nypanel, True, True)

    vrmin, vrmax = vr[0][0], vr[0][nvr]
    vtmin, vtmax = vt[0][0], vt[0][nvt]
    Emin, Emax = EE[0][0], EE[0][nE]
    Jmin, Jmax = JJ[0][0], JJ[0][nJ]

    # adjust marker size and etcetra
    # fs = fs_base / np.sqrt(nypanel)
    # tl = tl_base / np.sqrt(nypanel)
    # tw = tw_base / np.sqrt(nypanel)
    fs = fs_base
    tl = tl_base
    tw = tw_base


    for jj in range(nypanel):
        img_vel = ax_v[jj].imshow(vrvt[jj].T, extent = [vrmin, vrmax, vtmin, vtmax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = fvmin, vmax = fvmax), cmap = cmap_vel, aspect = 'auto', rasterized = True)
        img_ene = ax_E[jj].imshow(EJ[jj].T, extent = [Jmin, Jmax, Emin, Emax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = fEmin, vmax = fEmax), cmap = cmap_ene, aspect = 'auto', rasterized = True)



    col_grid = 'white'
    for ii in range(nxpanel * nypanel):
        ax_v[ii].set_xlim([vrmin, vrmax])
        ax_v[ii].set_ylim([vtmin, vtmax])
        ax_v[ii].tick_params(axis = 'both', direction = 'in', color = col_grid, bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
        ax_v[ii].spines['bottom'].set_color(col_grid)
        ax_v[ii].spines[   'top'].set_color(col_grid)
        ax_v[ii].spines[  'left'].set_color(col_grid)
        ax_v[ii].spines[ 'right'].set_color(col_grid)

        ax_E[ii].set_xlim([Jmin, Jmax])
        ax_E[ii].set_ylim([Emin, Emax])
        ax_E[ii].tick_params(axis = 'both', direction = 'in', color = col_grid, bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
        ax_E[ii].spines['bottom'].set_color(col_grid)
        ax_E[ii].spines[   'top'].set_color(col_grid)
        ax_E[ii].spines[  'left'].set_color(col_grid)
        ax_E[ii].spines[ 'right'].set_color(col_grid)
        ax_E[ii].yaxis.set_major_formatter(ticker.FuncFormatter(utils.scientific))

    col_frame = 'black'
    for jj in range(nypanel):
        ax_v[jj].spines['left'].set_color(col_frame)
        ax_v[jj].set_ylabel(r'$v_t$~(\si{km.s^{-1}})', fontsize = fs)
        ax_E[jj].spines['left'].set_color(col_frame)
        ax_E[jj].set_ylabel(r'$E$~(\si{erg.g^{-1}})', fontsize = fs)

        ax_v[(nxpanel - 1) * nypanel + jj].spines['right'].set_color(col_frame)
        ax_E[(nxpanel - 1) * nypanel + jj].spines['right'].set_color(col_frame)

    for ii in range(nxpanel):
        ax_v[ii * nypanel                ].spines['bottom'].set_color(col_frame)
        ax_v[ii * nypanel + (nypanel - 1)].spines[   'top'].set_color(col_frame)
        ax_v[ii * nypanel                ].set_xlabel(r'$v_r$~(\si{km.s^{-1}})', fontsize = fs)

        ax_E[ii * nypanel                ].spines['bottom'].set_color(col_frame)
        ax_E[ii * nypanel + (nypanel - 1)].spines[   'top'].set_color(col_frame)
        ax_E[ii * nypanel                ].set_xlabel(r'$J$~(\si{kpc.km.s^{-1}})', fontsize = fs)


    # add colorbar
    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax_v:
        xl, xr = at.get_position().x0, at.get_position().x1
        yb, yt = at.get_position().y0, at.get_position().y1

        if x0 > xl:
            x0 = xl
        if x1 < xr:
            x1 = xr
        if y0 > yb:
            y0 = yb
        if y1 < yt:
            y1 = yt
    cax_v = fig_v.add_axes([x1, y0, 0.05 / nxpanel, y1 - y0])
    cbar_v = fig_v.colorbar(img_vel, cax = cax_v, format=ticker.FuncFormatter(utils.scientific), label = r'$f$ (' + r'\si{M_\odot.km^{-2}.s^2}' + r')')
    cbar_v.solids.set_edgecolor('face')

    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax_E:
        xl, xr = at.get_position().x0, at.get_position().x1
        yb, yt = at.get_position().y0, at.get_position().y1

        if x0 > xl:
            x0 = xl
        if x1 < xr:
            x1 = xr
        if y0 > yb:
            y0 = yb
        if y1 < yt:
            y1 = yt
    cax_E = fig_E.add_axes([x1, y0, 0.05 / nxpanel, y1 - y0])
    cbar_E = fig_E.colorbar(img_ene, cax = cax_E, format=ticker.FuncFormatter(utils.scientific), label = r'$f$ (' + r'\si{M_\odot.kpc^{-1}.km^{-1}.s.erg^{-1}.g}' + r')')
    cbar_E.solids.set_edgecolor('face')


    # add caption
    col_caption = 'white'
    for ii in range(nxpanel):
        for jj in range(nypanel):
            caption  = '(' + '{:^c}'.format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ')'
            at_v = ax_v[ii * nypanel + jj]
            at_E = ax_E[ii * nypanel + jj]
            note = '~' + lab[nypanel - 1 - jj]

            at_v.text(0.03, 0.97, caption + note, color = col_caption, fontsize = fs, horizontalalignment = 'left', verticalalignment = 'top', transform = at_v.transAxes, bbox = dict(facecolor = get_cmap(cmap_vel)(0), edgecolor = 'None', alpha = 0.75))
            at_E.text(0.03, 0.97, caption + note, color = col_caption, fontsize = fs, horizontalalignment = 'left', verticalalignment = 'top', transform = at_E.transAxes, bbox = dict(facecolor = get_cmap(cmap_ene)(0), edgecolor = 'None', alpha = 0.75))


    # add current time
    if nypanel == 1:
        fig_v.suptitle(r'$t = {:.3f}$ {:<}'.format(time / 1000, 'Gyr'), x = 0.37, y = 1.0, fontsize = fs)
        fig_E.suptitle(r'$t = {:.3f}$ {:<}'.format(time / 1000, 'Gyr'), x = 0.37, y = 1.0, fontsize = fs)
    else:
        fig_v.suptitle(r'$t = {:.3f}$ {:<}'.format(time / 1000, 'Gyr'), y = 1.0, fontsize = fs)
        fig_E.suptitle(r'$t = {:.3f}$ {:<}'.format(time / 1000, 'Gyr'), y = 1.0, fontsize = fs)


    # save figures
    figname = 'fig/' + filename + '_m31phase'
    fig_v.savefig(figname + '_vel' + snapshot + '.png', format = 'png', dpi =  96, bbox_inches = 'tight')
    fig_E.savefig(figname + '_ene' + snapshot + '.png', format = 'png', dpi =  96, bbox_inches = 'tight')
    if outputPDF:
        fig_v.savefig(figname + '_vel' + snapshot + '.pdf', format = 'pdf', dpi = 300, bbox_inches = 'tight')
        fig_E.savefig(figname + '_ene' + snapshot + '.pdf', format = 'pdf', dpi = 300, bbox_inches = 'tight')
    plt.close('all')
