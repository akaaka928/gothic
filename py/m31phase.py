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


Nskip = 0
Nmbh = 0
init = 0
last = 264
last = 560
last = 428
last = 400
# Sigma_min, Sigma_max = 3.1e+4, 1.0e+8
# depth_min, depth_max = 3.1e+3, 3.1e+6
# phase_min, phase_max = 1.0e+4, 3.1e+6
lab = ['DM subhalo', 'star']



# preset_xmin, preset_xmax = -2.0, -8.0
# preset_ymin, preset_ymax = 1.0, 7.0
# preset_zmin, preset_zmax = m31.zm31 - 50.0, m31.zm31 + 200.0
# preset_vmin, preset_vmax = -480.0, -380.0


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


# cmap_Sigma = utils.generate_cmap(['darkblue', 'deepskyblue', 'lime', 'yellow', 'red', 'magenta', 'white'])
# cmap_depth = 'jet'
# cmap_phase = 'hot'
cmap_Sigma = 'plasma'
cmap_depth = 'cividis'
# cmap_depth = 'viridis'
cmap_phase = 'inferno'
# cmap_phase = 'viridis'


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
                EE[nypanel] = h5file[folder + 'EE']
                JJ[nypanel] = h5file[folder + 'JJ']

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
    fs = fs_base / np.sqrt(nypanel)
    tl = tl_base / np.sqrt(nypanel)
    tw = tw_base / np.sqrt(nypanel)


    for jj in range(nypanel):
        img_Sigma = ax[              jj].imshow(xy_map[jj].T, extent = [xmin, xmax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = Sigma_min, vmax = Sigma_max), cmap = cmap_Sigma, aspect = 'auto', rasterized = True)
        img_depth = ax[    nypanel + jj].imshow(yz_map[jj]  , extent = [zmin, zmax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = depth_min, vmax = depth_max), cmap = cmap_depth, aspect = 'auto', rasterized = True)
        img_phase = ax[2 * nypanel + jj].imshow(yv_map[jj]  , extent = [vmin, vmax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = phase_min, vmax = phase_max), cmap = cmap_phase, aspect = 'auto', rasterized = True)


    for jj in range(nypanel):
        # reference circles for projected distances
        c050 = plt.Circle((0, 0), m31.kpc2degree( 50.0), ec = col_ref, lw = lw_ref, fill = False)
        c100 = plt.Circle((0, 0), m31.kpc2degree(100.0), ec = col_ref, lw = lw_ref, fill = False)
        c150 = plt.Circle((0, 0), m31.kpc2degree(150.0), ec = col_ref, lw = lw_ref, fill = False)
        c200 = plt.Circle((0, 0), m31.kpc2degree(200.0), ec = col_ref, lw = lw_ref, fill = False)
        c250 = plt.Circle((0, 0), m31.kpc2degree(250.0), ec = col_ref, lw = lw_ref, fill = False)
        c300 = plt.Circle((0, 0), m31.kpc2degree(300.0), ec = col_ref, lw = lw_ref, fill = False)
        c350 = plt.Circle((0, 0), m31.kpc2degree(350.0), ec = col_ref, lw = lw_ref, fill = False)
        c400 = plt.Circle((0, 0), m31.kpc2degree(400.0), ec = col_ref, lw = lw_ref, fill = False)
        c450 = plt.Circle((0, 0), m31.kpc2degree(450.0), ec = col_ref, lw = lw_ref, fill = False)
        c500 = plt.Circle((0, 0), m31.kpc2degree(500.0), ec = col_ref, lw = lw_ref, fill = False)
        c550 = plt.Circle((0, 0), m31.kpc2degree(550.0), ec = col_ref, lw = lw_ref, fill = False)
        ax[jj].add_patch(c050)
        ax[jj].add_patch(c100)
        ax[jj].add_patch(c150)
        ax[jj].add_patch(c200)
        ax[jj].add_patch(c250)
        ax[jj].add_patch(c300)
        ax[jj].add_patch(c350)
        ax[jj].add_patch(c400)
        ax[jj].add_patch(c450)
        ax[jj].add_patch(c500)
        ax[jj].add_patch(c550)

        # reference ellipse of the M31 disk
        ax[              jj].plot(disk_xi[                    ::int((Ndisk - 1) / 2)], disk_eta[                    ::int((Ndisk - 1) / 2)], '-', color = col_M31disk, linewidth = lw_M31disk)# minor axis
        ax[              jj].plot(disk_xi[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_eta[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], '-', color = col_M31disk, linewidth = lw_M31disk)# major axis
        ax[              jj].plot(disk_xi                                            , disk_eta                                            , '-', color = col_M31disk, linewidth = lw_M31disk)
        ax[    nypanel + jj].plot(disk_D [                    ::int((Ndisk - 1) / 2)], disk_eta[                    ::int((Ndisk - 1) / 2)], '-', color = col_M31disk, linewidth = lw_M31disk)# minor axis
        ax[    nypanel + jj].plot(disk_D [int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_eta[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], '-', color = col_M31disk, linewidth = lw_M31disk)# major axis
        ax[    nypanel + jj].plot(disk_D                                             , disk_eta                                            , '-', color = col_M31disk, linewidth = lw_M31disk)

        # reference points of observation fields
        for kk in range(Nhsc):
            circle = plt.Circle((hsc_xi[kk], hsc_eta[kk]), hsc_rad[kk], ec = col_hsc, lw = lw_hsc, fill = False)
            ax[jj].add_patch(circle)

        # distance measurements to NW stream
        ax[              jj].plot(nws_s1_xi, nws_s1_eta, '-', color = col_nws, linewidth = lw_nws)
        ax[              jj].plot(nws_s2_xi, nws_s2_eta, '-', color = col_nws, linewidth = lw_nws)
        ax[              jj].plot(nws_s3_xi, nws_s3_eta, '-', color = col_nws, linewidth = lw_nws)
        ax[              jj].plot(nws_s4_xi, nws_s4_eta, '-', color = col_nws, linewidth = lw_nws)
        ax[    nypanel + jj].plot(nws_D , nws_eta, 's', color = col_nws, markerfacecolor = 'none', markersize = ms_nws)
        ax[    nypanel + jj].errorbar(nws_D, nws_eta, xerr = nws_Derr, ls = 'none', ecolor = col_nws, elinewidth = lw_nws)

        # line-of-sight velocity measurements of globular clusters
        ax[              jj].plot(nws_gc_xi , nws_gc_eta, 'o', color = col_gcs, markerfacecolor = 'none', markersize = ms_gcs)
        ax[2 * nypanel + jj].plot(nws_gc_vel, nws_gc_eta, 'o', color = col_gcs, markerfacecolor = 'none', markersize = ms_gcs)
        ax[2 * nypanel + jj].errorbar(nws_gc_vel, nws_gc_eta, xerr = nws_gc_verr, ls = 'none', ecolor = col_gcs, elinewidth = lw_gcs)


        # indicate wandering MBH location
        if Nmbh == 1:
            ax[              jj].plot(mbh_obs[0], mbh_obs[1], '+', color = col_wmbh, markersize = ms_wmbh)
            ax[    nypanel + jj].plot(mbh_obs[2], mbh_obs[1], '+', color = col_wmbh, markersize = ms_wmbh)


    if draw_whole_domain:
        draw_xmin, draw_xmax = xmin, xmax
        draw_ymin, draw_ymax = ymin, ymax
        draw_zmin, draw_zmax = zmin, zmax
        draw_vmin, draw_vmax = vmin, vmax
    else:
        draw_xmin, draw_xmax = preset_xmin, preset_xmax
        draw_ymin, draw_ymax = preset_ymin, preset_ymax
        draw_zmin, draw_zmax = preset_zmin, preset_zmax
        draw_vmin, draw_vmax = preset_vmin, preset_vmax


    for ii in range(nxpanel * nypanel):
        ax[ii].set_ylim([draw_ymin, draw_ymax])
        # ax[ii].set_yticks(ytics)
        ax[ii].tick_params(axis = 'both', direction = 'in', color = col_grid, bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
        ax[ii].spines['bottom'].set_color(col_grid)
        ax[ii].spines[   'top'].set_color(col_grid)
        ax[ii].spines[  'left'].set_color(col_grid)
        ax[ii].spines[ 'right'].set_color(col_grid)

    for jj in range(nypanel):
        ax[              jj].set_xlim([draw_xmin, draw_xmax])
        # ax[              jj].set_xticks(xtics)
        ax[              jj].spines['left'].set_color(col_frame)
        ax[              jj].set_ylabel(r'$\eta$ ({:<})'.format('degree'), fontsize = fs)
        ax[    nypanel + jj].set_xlim([draw_zmin, draw_zmax])
        # ax[    nypanel + jj].set_xticks(ztics)
        ax[2 * nypanel + jj].set_xlim([draw_vmin, draw_vmax])
        # ax[2 * nypanel + jj].set_xticks(vtics)
        ax[(nxpanel - 1) * nypanel + jj].spines['right'].set_color(col_frame)

    for ii in range(nxpanel):
        ax[ii * nypanel                ].spines['bottom'].set_color(col_frame)
        ax[ii * nypanel + (nypanel - 1)].spines[   'top'].set_color(col_frame)

    ax[0          ].set_xlabel(r'$\xi$ ({:<})'.format('degree'), fontsize = fs)
    ax[    nypanel].set_xlabel(r'$D$ ({:<})'.format('kpc'), fontsize = fs)
    ax[2 * nypanel].set_xlabel(r'$v_\mathrm{los}$ (km s$^{-1}$)', fontsize = fs)

    for ii in range(nxpanel):
        for jj in range(nypanel):
            caption  = '(' + '{:^c}'.format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ')'
            at = ax[ii * nypanel + jj]
            if ii == 0:
                cmap_here = cmap_Sigma
                note = ''
                if nypanel > 1:
                    note = '~' + lab[nypanel - 1 - jj]
            if ii == 1:
                cmap_here = cmap_depth
                note = ''
            if ii == 2:
                cmap_here = cmap_phase
                note = ''
            at.text(0.03, 0.97, caption + note, color = col_caption, fontsize = fs, horizontalalignment = 'left', verticalalignment = 'top', transform = at.transAxes, bbox = dict(facecolor = get_cmap(cmap_here)(0), edgecolor = 'None', alpha = 0.75))



    for ii in range(nxpanel):
        # add colorbar
        x0, x1 = 1, 0
        y0, y1 = 1, 0
        # for at in ax[ii * nypanel:ii * nypanel + (nypanel - 1)]:
        for at in ax[ii * nypanel:(ii + 1) * nypanel]:
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

        cax = fig.add_axes([x0, y1, x1 - x0, 0.05 / nypanel])
        if ii == 0:
            bar = fig.colorbar(img_Sigma, cax = cax, orientation = 'horizontal')
        if ii == 1:
            bar = fig.colorbar(img_depth, cax = cax, orientation = 'horizontal')
        if ii == 2:
            bar = fig.colorbar(img_phase, cax = cax, orientation = 'horizontal')
        bar.solids.set_edgecolor('face')
        if useDegree:
            if ii == 0:
                bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-2}})', fontsize = fs)
            if ii == 1:
                bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-1}.kpc^{-1}})', fontsize = fs)
            if ii == 2:
                bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-1}.km^{-1}.s})', fontsize = fs)
        else:
            if ii == 0:
                bar.set_label(r'$\Sigma$ (\si{M_\odot.kpc^{-2}})', fontsize = fs)
            if ii == 1:
                bar.set_label(r'$\Sigma$ (\si{M_\odot.kpc^{-2}})', fontsize = fs)
            if ii == 2:
                bar.set_label(r'$\Sigma$ (\si{M_\odot.kpc^{-1}.km^{-1}.s})', fontsize = fs)
        cax.tick_params(labelsize = fs)
        cax.xaxis.tick_top()
        cax.xaxis.set_label_position('top')


    # add current time
    if not monochrome:
        if nypanel == 1:
            fig.suptitle(r'$t = {:.3f}$ {:<}'.format(time / 1000, 'Gyr'), x = 0.37, y = 1.0, fontsize = fs)
        else:
            fig.suptitle(r'$t = {:.3f}$ {:<}'.format(time / 1000, 'Gyr'), y = 1.0, fontsize = fs)


    # save figures
    figname = 'fig/' + filename + '_map'
    if not draw_whole_domain:
        figname += '_hsc'
    if monochrome:
        figname += '_mono'
    figname += snapshot
    fig.savefig(figname + '.png', format = 'png', dpi =  96, bbox_inches = 'tight')
    if outputPDF:
        fig.savefig(figname + '.pdf', format = 'pdf', dpi = 300, bbox_inches = 'tight')
    plt.close('all')




import numpy as np
import math
import h5py

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm # for logarithmic plot in imshow
import matplotlib.ticker as ticker

import os.path as path
import multiprocessing as mp

import utils as utils


outputPDF = False

closeUp = False


fs_base = 24
tl_base = 6.0
tw_base = 1.0


senergy2astro = 1.0e+10


# filename = "m15ra1_4_run"
# init = 0
# last = 80
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["DM halo", "bulge"]
# fmin, fmax = 3.1e+0 / senergy2astro, 3.1e+3 / senergy2astro

# filename = "w3iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w3ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w3ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.5e+5 / senergy2astro

# filename = "w3ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 2.0e+5 / senergy2astro

# filename = "w3ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "w3ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+4 / senergy2astro, 1.0e+6 / senergy2astro

# filename = "w5iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w5ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w5ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w5ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 2.0e+5 / senergy2astro

# filename = "w5ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "w5ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 8.0e+5 / senergy2astro

# filename = "w7iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w7ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w7ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "w7ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "w7ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 3.1e+3 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "w7ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["GC"]
# fmin, fmax = 1.0e+4 / senergy2astro, 1.0e+6 / senergy2astro

# filename = "m12iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 4.0e+3 / senergy2astro

# filename = "m12ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 5.0e+3 / senergy2astro

filename = "m12ra1"
init = 0
last = 140
# last = 0
Nskip = 0
skip = [0]
lab = ["halo"]
fmin, fmax = 1.0e+2 / senergy2astro, 1.0e+4 / senergy2astro

# filename = "m12ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 1.0e+4 / senergy2astro

# filename = "m12ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 3.1e+4 / senergy2astro

# filename = "m12ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+0 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "m09iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "m09ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "m09ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "m09ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "m09ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "m09ra1_8"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 3.1e+5 / senergy2astro

# filename = "a00iso"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 3.1e+4 / senergy2astro

# filename = "a00ra2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 3.1e+4 / senergy2astro

# filename = "a00ra1"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 3.1e+4 / senergy2astro

# filename = "a00ra1_2"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 1.0e+2 / senergy2astro, 1.0e+5 / senergy2astro

# filename = "a00ra1_4"
# init = 0
# last = 140
# # last = 0
# Nskip = 0
# skip = [0]
# lab = ["halo"]
# fmin, fmax = 3.1e+2 / senergy2astro, 3.1e+5 / senergy2astro


Npt = 6
pt = ["o", "s", "^", "D", "x", "*"]
Nls = 4
ls = ["-", "--", ":", "-."]
Ncol, col = utils.set_color_palette_for_color_universal_design()


def draw_figure(fileid, kind, sphe):
    snapshot = "{:03d}".format(fileid)
    input_file = "dat/" + filename + ".plt" + snapshot + ".h5"

    # read snapshot
    h5file = h5py.File(input_file, "r")

    # read attributes
    time = h5file["/"].attrs["time"][0]

    nE = h5file["/"].attrs["nE"][0]
    nJ = h5file["/"].attrs["nJ"][0]

    # memory allocation for EJ-maps
    EJful = np.zeros((kind + 1, nJ, nE))
    EJsml = np.zeros((kind + 1, nJ, nE))
    EE    = np.zeros((kind + 1, nE + 1))
    Jful  = np.zeros((kind + 1, nJ + 1))
    Jsml  = np.zeros((kind + 1, nJ + 1))


    nxpanel = 1
    if closeUp:
        nxpanel = 2
    nypanel = 0
    Ndata = 0
    for ii in range(kind):
        num = h5file["attr" + str(ii)].attrs["number"][0]
        if (num > 1):
            if (Nskip == 0) or (ii not in skip):
                # EJ-maps
                folder = "EJplane" + str(ii) + "/"
                EJful[nypanel] = h5file[folder + "EJful"] / senergy2astro
                EJsml[nypanel] = h5file[folder + "EJsml"] / senergy2astro
                EE[nypanel] = h5file[folder + "E"] * senergy2astro
                Jful[nypanel] = h5file[folder + "Jful"]
                Jsml[nypanel] = h5file[folder + "Jsml"]

                nypanel += 1
                Ndata += 1

    # close the HDF5 file
    h5file.close()

    summary = False
    if nypanel > 1:
        for ii in range(nypanel):
            EJful[nypanel] += EJful[ii]
            EJsml[nypanel] += EJsml[ii]

        nypanel += 1
        summary = True


    EJful = np.maximum(EJful, fmin)
    EJful = np.minimum(EJful, fmax)
    EJsml = np.maximum(EJsml, fmin)
    EJsml = np.minimum(EJsml, fmax)


    # plot the data
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    Emin, Emax = EE[0][0], EE[0][nE]
    Jfulmin, Jfulmax = Jful[0][0], Jful[0][nJ]
    Jsmlmin, Jsmlmax = Jsml[0][0], Jsml[0][nJ]

    # adjust marker size and etcetra
    fs = fs_base / np.sqrt(max([nxpanel, nypanel]))
    tl = tl_base / np.sqrt(max([nxpanel, nypanel]))
    tw = tw_base / np.sqrt(max([nxpanel, nypanel]))

    if summary:
        img = ax[    nypanel - 1].imshow(EJful[nypanel - 1].T, extent = [Jfulmin, Jfulmax, Emin, Emax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        if closeUp:
            img = ax[2 * nypanel - 1].imshow(EJsml[nypanel - 1].T, extent = [Jsmlmin, Jsmlmax, Emin, Emax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

    for ii in range(Ndata):
        img = ax[ii          ].imshow(EJful[ii].T, extent = [Jfulmin, Jfulmax, Emin, Emax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")
        if closeUp:
            img = ax[ii + nypanel].imshow(EJsml[ii].T, extent = [Jsmlmin, Jsmlmax, Emin, Emax], origin = "lower", interpolation = "none", norm = LogNorm(vmin = fmin, vmax = fmax), cmap = my_cmap, aspect = "auto")

    for ii in range(nxpanel * nypanel):
        ax[ii].set_ylim([Emin, Emax])
        ax[ii].tick_params(axis = "both", direction = "in", color = "white", bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
        ax[ii].spines["bottom"].set_color("white")
        ax[ii].spines[   "top"].set_color("white")
        ax[ii].spines[  "left"].set_color("white")
        ax[ii].spines[ "right"].set_color("white")
        ax[ii].yaxis.set_major_formatter(ticker.FuncFormatter(utils.scientific))

    for ii in range(nypanel):
        ax[ii          ].set_xlim([Jfulmin, Jfulmax])
        ax[ii          ].spines[ "left"].set_color("black")
        ax[ii          ].set_ylabel(r"$E$~(\si{erg.g^{-1}})", fontsize = fs)
        if closeUp:
            ax[ii + nypanel].set_xlim([Jsmlmin, Jsmlmax])
            ax[ii + nypanel].spines["right"].set_color("black")

    for ii in range(nxpanel):
        ax[(ii + 1) * nypanel - 1].spines["top"].set_color("black")
        ax[ii * nypanel].spines["bottom"].set_color("black")
        ax[ii * nypanel].set_xlabel(r"$J$~(\si{kpc.km.s^{-1}})", fontsize = fs)


    # add colorbar
    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax:
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
    colorbar_ax = fig.add_axes([x1, y0, 0.05 / nxpanel, y1 - y0])
    cbar = fig.colorbar(img, cax = colorbar_ax, format=ticker.FuncFormatter(utils.scientific), label = r"$f$ (" + r"\si{M_\odot.kpc^{-1}.km^{-1}.s.erg^{-1}.g}" + r")")
    cbar.solids.set_edgecolor("face")

    # add current time
    fig.suptitle(r"$t = {:.1f}$".format(time / 1000) + r"~(\si{Gyr})", fontsize = fs)


    # save figures
    figname = "fig/" + filename + "_EJmap" + snapshot
    fig.savefig(figname + ".png", format = "png", dpi =  96, bbox_inches = "tight")
    if outputPDF:
        fig.savefig(figname + ".pdf", format = "pdf", dpi = 300, bbox_inches = "tight")
    plt.close("all")



def wrapper(argv):
    return draw_figure(*argv)


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rcParams['text.latex.preamble'] = [r"\usepackage{siunitx}"]

# set font size
plt.rcParams['font.size'] = 18
# plt.rcParams['font.size'] = 16

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# read number of all component(s)
txtfile = open("doc/" + filename + ".summary.txt", "r")
unit = int(txtfile.readline())
line = txtfile.readline()
item = line.split("\t")
Nkind = int(item[0])
Nsphe = int(item[1])
txtfile.close()


my_cmap = utils.generate_cmap(["darkblue", "deepskyblue", "lime", "yellow", "red", "magenta", "white"])
cores = int(np.ceil(mp.cpu_count() / 2))
pool = mp.Pool(cores)
args = [(ii, Nkind, Nsphe) for ii in range(init, last + 1, 1)]
pool.map(wrapper, args)
pool.close()
