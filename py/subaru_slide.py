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


smap_min, smap_max = 0.0, 0.0
vmap_min, vmap_max = 0.0, 0.0


from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


from argparse import ArgumentParser
def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-f', '--files', type=str,
                           help='Name of the target files')
    argparser.add_argument('-p', '--pdf',
                           action='store_true',
                           help='Whether to output figures in PDF format')
    argparser.add_argument('-m', '--monochrome',
                           action='store_true',
                           help='Whether to draw in monochrome')
    argparser.add_argument('-s', '--specify',
                           action='store_true',
                           help='Whether to specify the plot range')
    return argparser.parse_args()

# catch optional input arguments
args = get_option()
filename = args.files
outputPDF = args.pdf
monochrome = args.monochrome
specify = args.specify


fs_base = 16
tl_base = 6.0
tw_base = 1.0


init = 32
last = 112
# init = 40
# last = 40
# init = 56

# embed fonts
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

# use packages
plt.rcParams['text.latex.preamble'] = [r'\usepackage{physics,siunitx}']

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
else:
    kind = None
kind = comm.bcast(kind, root = 0)


cmap_map = 'plasma'
# cmap = 'cividis'
cmap_vel = 'inferno'
# cmap = 'viridis'


for fileid in range(init + rank, last + 1, size):
    snapshot = '{:03d}'.format(fileid)
    input_file = 'dat/' + filename + '.subaru' + snapshot + '.h5'

    # read snapshot
    with h5py.File(input_file, 'r') as h5file:
        # read attributes
        Nfield = h5file['subaru/'].attrs['Nfield'][0]
        Nmodel = h5file['subaru/'].attrs['Nmodel'][0] - 2
        Nnoise = h5file['subaru/'].attrs['Nnoise'][0]
        # Nfield = h5file['/'].attrs['Nfield'][0]
        # Nmodel = h5file['/'].attrs['Nmodel'][0]
        # Nnoise = h5file['/'].attrs['Nnoise'][0]

        field = h5file['subaru/'].attrs['field']
        model = h5file['subaru/'].attrs['model']
        SN = h5file['subaru/'].attrs['SN']
        # field = h5file['/'].attrs['field']
        # model = h5file['/'].attrs['model']

        nx = h5file['/'].attrs['nx'][0]
        ny = h5file['/'].attrs['ny'][0]
        nv = h5file['/'].attrs['nv'][0]

        # memory allocation for surface density maps
        xy_noise = np.zeros((kind, Nfield, Nnoise, nx, ny))
        xy_map = np.zeros((kind, Nfield, Nmodel, nx, ny))
        xx     = np.zeros((kind, Nfield, nx + 1))
        yy     = np.zeros((kind, Nfield, ny + 1))

        # memory allocation for phase-space densty map
        yv_noise = np.zeros((kind, Nfield, Nnoise, ny, nv))
        yv_map = np.zeros((kind, Nfield, Nmodel, ny, nv))
        vv     = np.zeros((kind, Nfield, nv + 1))

        for kk in range(kind):
            folder = 'component' + str(kk) + '/'
            for ff in range(Nfield):
                # print(folder)
                # print(field[ff].decode('utf-8'))
                xx[kk][ff] = h5file[folder + field[ff].decode('utf-8') + '-xi']
                yy[kk][ff] = h5file[folder + field[ff].decode('utf-8') + '-eta']
                vv[kk][ff] = h5file[folder + field[ff].decode('utf-8') + '-vlos']
                xy_map[kk][ff][0] = h5file[folder + field[ff].decode('utf-8') + '-' + model[0].decode('utf-8') + '-map']
                yv_map[kk][ff][0] = h5file[folder + field[ff].decode('utf-8') + '-' + model[0].decode('utf-8') + '-vel']
                for mm in range(Nmodel - 1):
                    xy_map[kk][ff][1 + mm] = h5file[folder + field[ff].decode('utf-8') + '-' + model[1 + 2 + mm].decode('utf-8') + '-map']
                    yv_map[kk][ff][1 + mm] = h5file[folder + field[ff].decode('utf-8') + '-' + model[1 + 2 + mm].decode('utf-8') + '-vel']
                for mm in range(Nnoise):
                    xy_noise[kk][ff][mm] = h5file[folder + field[ff].decode('utf-8') + '-SN' + '{:02d}'.format(SN[mm]) + '-map']
                    yv_noise[kk][ff][mm] = h5file[folder + field[ff].decode('utf-8') + '-SN' + '{:02d}'.format(SN[mm]) + '-vel']
        # close the HDF5 file


    # xy_map = np.maximum(xy_map, Sigma_min)
    # xy_map = np.minimum(xy_map, Sigma_max)
    # yv_map = np.maximum(yv_map, phase_min)
    # yv_map = np.minimum(yv_map, phase_max)

    nxpanel = Nmodel
    nypanel = 2 # w/o or w/z noise
    smap = np.zeros((nxpanel, nypanel, nx, ny))
    vmap = np.zeros((nxpanel, nypanel, nv, ny))

    for kk in range(kind):
        for ff in range(Nfield):

            for ii in range(Nmodel):
                smap[ii][nypanel - 1] = xy_map[kk][ff][ii]
                vmap[ii][nypanel - 1] = yv_map[kk][ff][ii]

            # add noise component
            for noise in range(Nnoise):
                for ii in range(Nmodel):
                    smap[ii][0] = xy_map[kk][ff][ii] + xy_noise[kk][ff][noise]
                    vmap[ii][0] = yv_map[kk][ff][ii] + yv_noise[kk][ff][noise]
                    # smap[ii][0] = xy_noise[kk][ff][noise]
                    # vmap[ii][0] = yv_noise[kk][ff][noise]

                if not specify:
                    smap_max = np.amax(smap)
                    # tmp = np.ma.masked_equal(smap, 0.0, copy = False)
                    # smap_min = np.amin(tmp)
                    smap_min = smap_max * 1.0e-2
                    vmap_max = np.amax(vmap)
                    # tmp = np.ma.masked_equal(vmap, 0.0, copy = False)
                    # vmap_min = np.amin(tmp)
                    vmap_min = vmap_max * 1.0e-2

                smap = np.maximum(smap, smap_min)
                smap = np.minimum(smap, smap_max)
                vmap = np.maximum(vmap, vmap_min)
                vmap = np.minimum(vmap, vmap_max)

                xmin, xmax = xx[kk][ff][0], xx[kk][ff][nx]
                ymin, ymax = yy[kk][ff][0], yy[kk][ff][ny]
                vmin, vmax = vv[kk][ff][0], vv[kk][ff][nv]

                # adjust marker size and etcetra
                fs = fs_base# / np.sqrt(nypanel)
                tl = tl_base# / np.sqrt(nypanel)
                tw = tw_base# / np.sqrt(nypanel)

                fig_map = utils.set_figure(nxpanel, nypanel)
                fig_vel = utils.set_figure(nxpanel, nypanel)
                ax_map = [0] * nxpanel * nypanel
                ax_vel = [0] * nxpanel * nypanel
                utils.locate_panels(fig_map, ax_map, nxpanel, nypanel, True, True)
                utils.locate_panels(fig_vel, ax_vel, nxpanel, nypanel, True, True)

                for ii in range(nxpanel):
                    for jj in range(nypanel):
                        img_map = ax_map[jj + ii * nypanel].imshow(smap[ii][jj].T, extent = [xmin, xmax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = smap_min, vmax = smap_max), cmap = cmap_map, aspect = 'auto', rasterized = True)
                        img_vel = ax_vel[jj + ii * nypanel].imshow(vmap[ii][jj]  , extent = [vmin, vmax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = vmap_min, vmax = vmap_max), cmap = cmap_vel, aspect = 'auto', rasterized = True)


                for at in ax_map:
                    at.set_xlim([xmax, xmin])
                    at.set_ylim([ymin, ymax])
                    at.tick_params(axis = 'both', direction = 'in', color = 'white', bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
                    at.spines['bottom'].set_color('white')
                    at.spines[   'top'].set_color('white')
                    at.spines[  'left'].set_color('white')
                    at.spines[ 'right'].set_color('white')
                for at in ax_vel:
                    at.set_xlim([vmin, vmax])
                    at.set_ylim([ymin, ymax])
                    at.tick_params(axis = 'both', direction = 'in', color = 'white', bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
                    at.spines['bottom'].set_color('white')
                    at.spines[   'top'].set_color('white')
                    at.spines[  'left'].set_color('white')
                    at.spines[ 'right'].set_color('white')

                for ii in range(nxpanel):
                    ax_map[                ii * nypanel].spines['bottom'].set_color('black')
                    ax_map[(nypanel - 1) + ii * nypanel].spines[   'top'].set_color('black')
                    ax_vel[                ii * nypanel].spines['bottom'].set_color('black')
                    ax_vel[(nypanel - 1) + ii * nypanel].spines[   'top'].set_color('black')

                    ax_map[                ii * nypanel].set_xlabel(r'$\xi$~(\si{deg})', fontsize = fs)
                    ax_vel[                ii * nypanel].set_xlabel(r'$v_\mathrm{los}$~(\si{km.s^{-1}})', fontsize = fs)

                for jj in range(nypanel):
                    ax_map[jj                          ].spines[ 'left'].set_color('black')
                    ax_map[jj + (nxpanel - 1) * nypanel].spines['right'].set_color('black')
                    ax_vel[jj                          ].spines[ 'left'].set_color('black')
                    ax_vel[jj + (nxpanel - 1) * nypanel].spines['right'].set_color('black')

                    ax_map[jj                          ].set_ylabel(r'$\eta$~(\si{deg})', fontsize = fs)
                    ax_vel[jj                          ].set_ylabel(r'$\eta$~(\si{deg})', fontsize = fs)

                for ii in range(nxpanel):
                    for jj in range(nypanel):
                        caption = '(' + '{:^c}'.format(97 + ii + nxpanel * (nypanel - 1 - jj)) + ')'
                        if ii == 0:
                            caption += '~w/o encounter'
                        if ii == 1:
                            caption += r'~$M = \SI{e+8}{M_\odot}$'
                        if ii == 2:
                            caption += r'~$M = \SI{e+8.5}{M_\odot}$'
                        if ii == 3:
                            caption += r'~$M = \SI{e+9}{M_\odot}$'
                        if ii == 4:
                            caption += r'~$M = \SI{e+9.5}{M_\odot}$'
                        if jj == 0:
                            caption += r'~$+$~noise'
                        at = ax_map[jj + ii * nypanel]
                        at.text(0.03, 0.97, caption, color = 'white', fontsize = fs, horizontalalignment = 'left', verticalalignment = 'top', transform = at.transAxes, bbox = dict(facecolor = get_cmap(cmap_map)(0), edgecolor = 'None', alpha = 0.75))
                        at = ax_vel[jj + ii * nypanel]
                        at.text(0.03, 0.97, caption, color = 'white', fontsize = fs, horizontalalignment = 'left', verticalalignment = 'top', transform = at.transAxes, bbox = dict(facecolor = get_cmap(cmap_vel)(0), edgecolor = 'None', alpha = 0.75))


                # add colorbar
                x0, x1 = 1, 0
                y0, y1 = 1, 0
                for at in ax_map:
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
                cax = fig_map.add_axes([x1, y0, 0.05 / nxpanel, y1 - y0])
                bar = fig_map.colorbar(img_map, cax = cax)
                bar.solids.set_edgecolor('face')
                bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-2}})', fontsize = fs)
                cax.tick_params(labelsize = fs)

                x0, x1 = 1, 0
                y0, y1 = 1, 0
                for at in ax_vel:
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
                cax = fig_vel.add_axes([x1, y0, 0.05 / nxpanel, y1 - y0])
                bar = fig_vel.colorbar(img_vel, cax = cax)
                bar.solids.set_edgecolor('face')
                bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-1}.km^{-1}.s})', fontsize = fs)
                cax.tick_params(labelsize = fs)

                # save figures
                pos = model[Nmodel - 1].decode('utf-8').rfind('-')
                figname = 'nws' + model[Nmodel - 1].decode('utf-8')[pos:] + '-kind' + str(kk) + '-' + field[ff].decode('utf-8') + '-SN' + '{:02d}'.format(SN[noise])
                fig_map.savefig('fig/' + figname + '-map' + snapshot + '.png', format = 'png', dpi = 96, bbox_inches = 'tight')
                fig_vel.savefig('fig/' + figname + '-vel' + snapshot + '.png', format = 'png', dpi = 96, bbox_inches = 'tight')
                if outputPDF:
                    fig_map.savefig('fig/' + figname + '-map' + snapshot + '.pdf', format = 'pdf', dpi = 300, bbox_inches = 'tight')
                    fig_vel.savefig('fig/' + figname + '-vel' + snapshot + '.pdf', format = 'pdf', dpi = 300, bbox_inches = 'tight')


                plt.close('all')
