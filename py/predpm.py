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
import m31 as m31


from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


from argparse import ArgumentParser
def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-f', '--files', type=str,
                           help='Name of the target files')
    argparser.add_argument('-r', '--ref', type=str,
                           default='../cont/dat/nws-continue',
                           help='Name of the target files')
    argparser.add_argument('-c', '--continued',
                           action='store_true',
                           help='Visualize continued simulations')
    argparser.add_argument('-s', '--specify',
                           action='store_false',
                           help='Whether to set the plot domain to cover whole evolution of the system (e.g., merger simulations)')
    argparser.add_argument('-p', '--pdf',
                           action='store_true',
                           help='Whether to output figures in PDF format')
    return argparser.parse_args()

# catch optional input arguments
args = get_option()
filename = args.files
refname = args.ref
continued = args.continued
draw_whole_domain = args.specify
outputPDF = args.pdf



fs_base = 16
tl_base = 6.0
tw_base = 1.0

col_ref = '#84919e'
lw_ref_base = 1.0

col_M31disk = 'white'
lw_M31disk_base = 1.0


pt = ['o', 's', '^', 'D']
ls = ['-', ':', '-.', '--']
Ncol, col = utils.set_color_palette_for_color_universal_design()

col_hsc, col_nws, col_gcs = 'white', col[4], col[5]
lw_hsc_base, lw_nws_base, lw_gcs_base = 1.5, 1.5, 1.5
ms_hsc_base, ms_nws_base, ms_gcs_base = 1, 4, 4

col_frame = 'black'
col_grid  = 'white'
col_caption = 'white'


time_offset = 0.0
if continued:
    # time_offset = -775.0
    # time_offset = -4875.0
    time_offset = -1000.0
    # time_offset = -1200.0


Nskip = 0
init = 0
# last = 264
# last = 428
last = 560
Sigma_min, Sigma_max = 3.1e+4, 1.0e+8
phase_min, phase_max = 1.0e+4, 3.1e+6
lab = ['DM subhalo', 'star']



preset_xmin, preset_xmax = -2.0, -8.0
preset_ymin, preset_ymax = 1.0, 7.0
preset_vxmin, preset_vxmax = -480.0, -380.0
preset_vymin, preset_vymax = -480.0, -380.0
# xtics = [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0]
# ytics = [-2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
# ztics = [750.0, 800.0, 850.0, 900.0, 950.0]


# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# use packages
plt.rcParams['text.latex.preamble'] = r'\usepackage{physics,siunitx}'

# set font size
plt.rcParams['font.size'] = fs_base

# specify direction of ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# set reference ellipse of M31 disk
Ndisk, disk_xi, disk_eta, disk_D = m31.disk_ellipse()

# set reference points of observations by Komiyama et al. (2018)
Nhsc, hsc_xi, hsc_eta, hsc_rad, Nnws, nws_xi, nws_eta, nws_D, nws_Dep, nws_Dem, nws_s1_xi, nws_s1_eta, nws_s2_xi, nws_s2_eta, nws_s3_xi, nws_s3_eta, nws_s4_xi, nws_s4_eta = m31.NWS_Komiyama2018()
nws_Derr = np.zeros((2, Nnws))
for ii in range(Nnws):
    nws_Derr[0][ii] = nws_Dem[ii]
    nws_Derr[1][ii] = nws_Dep[ii]

# set reference points of line-of-sight velocity measurements of GCs by Veljanoski et al. (2014)
Ngcs, nws_gc_xi, nws_gc_eta, nws_gc_vel, nws_gc_vep, nws_gc_vem = m31.NWS_velocity()
nws_gc_verr = np.zeros((2, Ngcs))
for ii in range(Ngcs):
    nws_gc_verr[0][ii] = nws_gc_vem[ii]
    nws_gc_verr[1][ii] = nws_gc_vep[ii]

cmap_Sigma = 'plasma'
cmap_phase = 'inferno'


for fileid in range(init + rank, last + 1, size):
    snapshot = '{:03d}'.format(fileid)
    obs_ref = refname + '.m31obs' + snapshot + '.h5'
    prd_ref = refname + '.m31prd' + snapshot + '.h5'
    obs_dat = 'dat/' + filename + '.m31obs' + snapshot + '.h5'
    prd_dat = 'dat/' + filename + '.m31prd' + snapshot + '.h5'

    nxpanel = 3

    # read reference snapshot
    with h5py.File(obs_ref, 'r') as h5file:
        # read attributes
        time = h5file['/'].attrs['time'][0] + time_offset
        useDegree = h5file['/'].attrs['useDegree'][0]

        kind = h5file['/'].attrs['kinds'][0]
        nx = h5file['/'].attrs['nx'][0]
        ny = h5file['/'].attrs['ny'][0]

        # memory allocation for surface density maps
        xy_map = np.zeros((2 * kind, nx, ny))
        xx     = np.zeros((2 * kind, nx + 1))
        yy     = np.zeros((2 * kind, ny + 1))

        nypanel = 0
        Ndata = 0
        for ii in range(kind):
            if (Nskip == 0) or (ii not in skip):
                # read surface density maps
                folder = 'field' + str(ii) + '/'
                xy_map[2 * kind - 1 - nypanel] = h5file[folder + 'Sigma_xy']
                xx[2 * kind - 1 - nypanel] = h5file[folder +   'xi']
                yy[2 * kind - 1 - nypanel] = h5file[folder +  'eta']

                nypanel += 1
                Ndata += 1
        # close the HDF5 file

    # read snapshot
    with h5py.File(obs_dat, 'r') as h5file:
        for ii in range(kind):
            if (Nskip == 0) or (ii not in skip):
                # read surface density maps
                folder = 'field' + str(ii) + '/'
                xy_map[2 * kind - 1 - nypanel] = h5file[folder + 'Sigma_xy']
                xx[2 * kind - 1 - nypanel] = h5file[folder +   'xi']
                yy[2 * kind - 1 - nypanel] = h5file[folder +  'eta']

                nypanel += 1
                Ndata += 1
        # close the HDF5 file




    # read reference snapshot
    with h5py.File(prd_ref, 'r') as h5file:
        # read attributes
        time = h5file['/'].attrs['time'][0] + time_offset
        useDegree = h5file['/'].attrs['useDegree'][0]

        kind = h5file['/'].attrs['kinds'][0]
        ny = h5file['/'].attrs['ny'][0]
        nvx = h5file['/'].attrs['nvxi'][0]
        nvy = h5file['/'].attrs['nveta'][0]

        # memory allocation for surface density maps
        yy_vx = np.zeros((2 * kind, ny, nvx))
        yy_vy = np.zeros((2 * kind, ny, nvy))
        eta   = np.zeros((2 * kind, ny + 1))
        vx    = np.zeros((2 * kind, nvx + 1))
        vy    = np.zeros((2 * kind, nvy + 1))

        nypanel = 0
        Ndata = 0
        for ii in range(kind):
            if (Nskip == 0) or (ii not in skip):
                # read surface density maps
                folder = 'field' + str(ii) + '/'
                yy_vx[2 * kind - 1 - nypanel] = h5file[folder + 'eta-vxi']
                yy_vy[2 * kind - 1 - nypanel] = h5file[folder + 'eta-veta']
                eta[2 * kind - 1 - nypanel] = h5file[folder +  'eta']
                vx[2 * kind - 1 - nypanel] = h5file[folder +  'vxi']
                vy[2 * kind - 1 - nypanel] = h5file[folder +  'veta']

                nypanel += 1
                Ndata += 1
        # close the HDF5 file

    # read snapshot
    with h5py.File(prd_dat, 'r') as h5file:
        for ii in range(kind):
            if (Nskip == 0) or (ii not in skip):
                # read surface density maps
                folder = 'field' + str(ii) + '/'
                yy_vx[2 * kind - 1 - nypanel] = h5file[folder + 'eta-vxi']
                yy_vy[2 * kind - 1 - nypanel] = h5file[folder + 'eta-veta']
                eta[2 * kind - 1 - nypanel] = h5file[folder +  'eta']
                vx[2 * kind - 1 - nypanel] = h5file[folder +  'vxi']
                vy[2 * kind - 1 - nypanel] = h5file[folder +  'veta']

                nypanel += 1
                Ndata += 1
        # close the HDF5 file


    xy_map = np.maximum(xy_map, Sigma_min)
    xy_map = np.minimum(xy_map, Sigma_max)
    yy_vx  = np.maximum(yy_vx , phase_min)
    yy_vx  = np.minimum(yy_vx , phase_max)
    yy_vy  = np.maximum(yy_vy , phase_min)
    yy_vy  = np.minimum(yy_vy , phase_max)


    # plot the data
    fig = utils.set_figure(nxpanel, nypanel)
    ax = [0] * nxpanel * nypanel
    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

    xmin, xmax = xx[0][0], xx[0][nx]
    ymin, ymax = yy[0][0], yy[0][ny]
    vxmin, vxmax = vx[0][0], vx[0][nvx]
    vymin, vymax = vy[0][0], vy[0][nvy]


    # adjust marker size and etcetra
    fs = fs_base / np.sqrt(nypanel)
    tl = tl_base / np.sqrt(nypanel)
    tw = tw_base / np.sqrt(nypanel)
    lw_ref = lw_ref_base / nypanel
    lw_M31disk = lw_M31disk_base / nypanel
    lw_hsc = lw_hsc_base / nypanel
    lw_nws = lw_nws_base / nypanel
    ms_hsc = ms_hsc_base / nypanel
    ms_nws = ms_nws_base / nypanel



    for jj in range(nypanel):
        img_xxyy = ax[              jj].imshow(xy_map[jj].T, extent = [xmin, xmax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = Sigma_min, vmax = Sigma_max), cmap = cmap_Sigma, aspect = 'auto', rasterized = True)
        img_yyvx = ax[    nypanel + jj].imshow(yy_vx[jj], extent = [vxmin, vymax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = phase_min, vmax = phase_max), cmap = cmap_phase, aspect = 'auto', rasterized = True)
        img_yyvy = ax[2 * nypanel + jj].imshow(yy_vy[jj], extent = [vymin, vymax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = phase_min, vmax = phase_max), cmap = cmap_phase, aspect = 'auto', rasterized = True)


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
        ax[jj].plot(disk_xi[                    ::int((Ndisk - 1) / 2)], disk_eta[                    ::int((Ndisk - 1) / 2)], '-', color = col_M31disk, linewidth = lw_M31disk)# minor axis
        ax[jj].plot(disk_xi[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], disk_eta[int((Ndisk - 1) / 4)::int((Ndisk - 1) / 2)], '-', color = col_M31disk, linewidth = lw_M31disk)# major axis
        ax[jj].plot(disk_xi                                            , disk_eta                                            , '-', color = col_M31disk, linewidth = lw_M31disk)

        # reference points of observation fields
        for kk in range(Nhsc):
            circle = plt.Circle((hsc_xi[kk], hsc_eta[kk]), hsc_rad[kk], ec = col_hsc, lw = lw_hsc, fill = False)
            ax[jj].add_patch(circle)

        # distance measurements to NW stream
        ax[jj].plot(nws_s1_xi, nws_s1_eta, '-', color = col_nws, linewidth = lw_nws)
        ax[jj].plot(nws_s2_xi, nws_s2_eta, '-', color = col_nws, linewidth = lw_nws)
        ax[jj].plot(nws_s3_xi, nws_s3_eta, '-', color = col_nws, linewidth = lw_nws)
        ax[jj].plot(nws_s4_xi, nws_s4_eta, '-', color = col_nws, linewidth = lw_nws)


    if draw_whole_domain:
        draw_xmin, draw_xmax = xmin, xmax
        draw_ymin, draw_ymax = ymin, ymax
        draw_vxmin, draw_vxmax = vxmin, vxmax
        draw_vymin, draw_vymax = vymin, vymax
    else:
        draw_xmin, draw_xmax = preset_xmin, preset_xmax
        draw_ymin, draw_ymax = preset_ymin, preset_ymax
        draw_vxmin, draw_vxmax = preset_vxmin, preset_vxmax
        draw_vymin, draw_vymax = preset_vymin, preset_vymax


    for ii in range(nxpanel * nypanel):
        ax[ii].set_ylim([draw_ymin, draw_ymax])
        ax[ii].tick_params(axis = 'both', direction = 'in', color = col_grid, bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
        ax[ii].spines['bottom'].set_color(col_grid)
        ax[ii].spines[   'top'].set_color(col_grid)
        ax[ii].spines[  'left'].set_color(col_grid)
        ax[ii].spines[ 'right'].set_color(col_grid)

    for jj in range(nypanel):
        ax[              jj].set_xlim([draw_xmin, draw_xmax])
        ax[              jj].spines['left'].set_color(col_frame)
        ax[              jj].set_ylabel(r'$\eta$~($\si{deg}$)', fontsize = fs)
        ax[    nypanel + jj].set_xlim([draw_vxmin, draw_vxmax])
        ax[2 * nypanel + jj].set_xlim([draw_vymin, draw_vymax])
        ax[(nxpanel - 1) * nypanel + jj].spines['right'].set_color(col_frame)

    for ii in range(nxpanel):
        ax[ii * nypanel                ].spines['bottom'].set_color(col_frame)
        ax[ii * nypanel + (nypanel - 1)].spines[   'top'].set_color(col_frame)

    ax[0          ].set_xlabel(r'$\xi$~($\si{deg}$)', fontsize = fs)
    ax[    nypanel].set_xlabel(r'$v_\xi$~($\si{km.s^{-1}}$)', fontsize = fs)
    ax[2 * nypanel].set_xlabel(r'$v_\eta$~($\si{km.s^{-1}}$)', fontsize = fs)

    for ii in range(nxpanel):
        for jj in range(nypanel):
            caption  = '(' + '{:^c}'.format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ')'
            at = ax[ii * nypanel + jj]
            if ii == 0:
                cmap_here = cmap_Sigma
                note = ''
                if kind > 1:
                    note = '~' + lab[kind - 1 - (jj % kind)]
                if jj < kind:
                    note += '~' + 'w/ DM subhalo perturbation'
                else:
                    note += '~' + 'w/o DM subhalo perturbation'
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
                bar.set_label(r'$\Sigma$~(\si{M_\odot.deg^{-2}})', fontsize = fs)
            else:
                bar.set_label(r'$\Sigma$~(\si{M_\odot.deg^{-1}.km^{-1}.s})', fontsize = fs)
        else:
            if ii == 0:
                bar.set_label(r'$\Sigma$~(\si{M_\odot.kpc^{-2}})', fontsize = fs)
            else:
                bar.set_label(r'$\Sigma$~(\si{M_\odot.kpc^{-1}.km^{-1}.s})', fontsize = fs)
        cax.tick_params(labelsize = fs)
        cax.xaxis.tick_top()
        cax.xaxis.set_label_position('top')


    # add current time
    if nypanel == 1:
        fig.suptitle(r'$t = {:.3f}$ {:<}'.format(time / 1000, 'Gyr'), x = 0.37, y = 1.0, fontsize = fs)
    else:
        fig.suptitle(r'$t = {:.3f}$ {:<}'.format(time / 1000, 'Gyr'), y = 1.0, fontsize = fs)


    # save figures
    figname = 'fig/' + filename + '_pm'
    if not draw_whole_domain:
        figname += '_hsc'
    figname += snapshot
    fig.savefig(figname + '.png', format = 'png', dpi =  96, bbox_inches = 'tight')
    if outputPDF:
        fig.savefig(figname + '.pdf', format = 'pdf', dpi = 300, bbox_inches = 'tight')
    plt.close('all')
