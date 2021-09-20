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


from argparse import ArgumentParser
def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-f', '--files', type=str,
                           help='Name of the target files')
    argparser.add_argument('-p', '--pdf',
                           action='store_true',
                           help='Whether to output figures in PDF format')
    return argparser.parse_args()

# catch optional input arguments
args = get_option()
filename = args.files
outputPDF = args.pdf


fs_base = 16
tl_base = 6.0
tw_base = 1.0

col_ref = '#84919e'
lw_ref_base = 1.0

col_M31disk = 'white'
lw_M31disk_base = 1.0

col_wmbh = 'black'
ms_wmbh_base = 3


pt = ['o', 's', '^', 'D']
ls = ['-', ':', '-.', '--']
Ncol, col = utils.set_color_palette_for_color_universal_design()

col_hsc, col_nws, col_gcs = 'white', col[4], col[5]
lw_hsc_base, lw_nws_base, lw_gcs_base = 1.5, 1.5, 1.5
ms_hsc_base, ms_nws_base, ms_gcs_base = 1, 4, 4

col_frame = 'black'
col_grid  = 'white'
col_caption = 'white'


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


view_cmap = 'plasma'
dist_cmap = 'cividis'
# cmap_depth = 'viridis'
# cmap_phase = 'inferno'




with h5py.File('dat/' + filename + '.h5', 'r') as h5file:
    # read attributes
    Nx = h5file['data/'].attrs['MAP_NX'][0]
    Ny = h5file['data/'].attrs['MAP_NY'][0]
    # Nf = h5file['data/'].attrs['DISTANCE_NFIELD'][0]
    # Nd = h5file['data/'].attrs['DISTANCE_NDEPTH'][0]
    best = h5file['info/'].attrs['best_index'][0]
    Xmin = h5file['NWstream/'].attrs['Xmin'][0]
    Xmax = h5file['NWstream/'].attrs['Xmax'][0]
    Ymin = h5file['NWstream/'].attrs['Ymin'][0]
    Ymax = h5file['NWstream/'].attrs['Ymax'][0]
    axis = h5file['NWstream/'].attrs['axis'][0]
    fwhm = h5file['NWstream/'].attrs['width'][0]
    gYmin = h5file['NWstream/'].attrs['gradient_ymin'][0]
    gYmax = h5file['NWstream/'].attrs['gradient_ymax'][0]
    # print(axis)
    # print(fwhm)

    # # memory allocation
    # view = np.zeros((Ny, Nx))
    # dist = np.zeros((Nf, Nd))

    folder = 'data/phi_' + str(best) + '/'
    view = np.array(h5file[folder + 'map'])
    dist = np.array(h5file[folder + 'box'])

    # close the HDF5 file

mass2astro = 1.0e+8
dX, dY = (Xmax - Xmin) / float(Nx), (Ymax - Ymin) / float(Ny)
dSinv = mass2astro / (dX * dY)
view *= dSinv
dist *= mass2astro

view_max = np.amax(view)
view_tmp = np.ma.masked_equal(view, 0.0, copy = False)
view_min = np.amin(view_tmp)
dist_max = np.amax(dist)
dist_tmp = np.ma.masked_equal(dist, 0.0, copy = False)
dist_min = np.amin(dist_tmp)

Ndat = 2
axis_y = np.linspace(Ymin, Ymax, Ndat)
axis_x = np.linspace(axis, axis, Ndat)
edge_l = np.linspace(axis - 0.5 * fwhm, axis - 0.5 * fwhm, Ndat)
edge_r = np.linspace(axis + 0.5 * fwhm, axis + 0.5 * fwhm, Ndat)
edge_L = np.linspace(axis -       fwhm, axis -       fwhm, Ndat)
edge_R = np.linspace(axis +       fwhm, axis +       fwhm, Ndat)
area_x = np.linspace(Xmin, Xmax, Ndat)
area_b = np.linspace(gYmin, gYmin, Ndat)
area_t = np.linspace(gYmax, gYmax, Ndat)


# plot the data
nxpanel, nypanel = 1, 1
map_fig = utils.set_figure(nxpanel, nypanel)
map_ax = [0] * nxpanel * nypanel
utils.locate_panels(map_fig, map_ax, nxpanel, nypanel, True, True)
box_fig = utils.set_figure(nxpanel, nypanel)
box_ax = [0] * nxpanel * nypanel
utils.locate_panels(box_fig, box_ax, nxpanel, nypanel, True, True)



# adjust marker size and etcetra
fs = fs_base / np.sqrt(nypanel)
tl = tl_base / np.sqrt(nypanel)
tw = tw_base / np.sqrt(nypanel)
lw_ref = lw_ref_base / nypanel


map_img = map_ax[0].imshow(view, extent = [Xmin, Xmax, Ymin, Ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = view_min, vmax = view_max), cmap = view_cmap, aspect = 'equal', rasterized = True)
# map_img = map_ax[0].imshow(view, extent = [Xmin, Xmax, Ymin, Ymax], origin = 'lower', interpolation = 'none', cmap = view_cmap, aspect = 'equal', rasterized = True)
box_img = box_ax[0].imshow(dist, origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = dist_min, vmax = dist_max), cmap = dist_cmap, aspect = 'equal', rasterized = True)
# box_img = box_ax[0].imshow(dist, origin = 'lower', interpolation = 'none', cmap = dist_cmap, aspect = 'equal', rasterized = True)

map_ax[0].plot(axis_x, axis_y, linestyle = '-', color = 'black')
map_ax[0].plot(edge_l, axis_y, linestyle = '--', color = 'black')
map_ax[0].plot(edge_r, axis_y, linestyle = '--', color = 'black')
map_ax[0].plot(edge_L, axis_y, linestyle = ':', color = 'black')
map_ax[0].plot(edge_R, axis_y, linestyle = ':', color = 'black')
map_ax[0].plot(area_x, area_b, linestyle = '-.', color = 'black')
map_ax[0].plot(area_x, area_t, linestyle = '-.', color = 'black')


for at in map_ax:
    at.tick_params(axis = 'both', direction = 'in', color = col_grid, bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
    # at.spines['bottom'].set_color(col_grid)
    # at.spines[   'top'].set_color(col_grid)
    # at.spines[  'left'].set_color(col_grid)
    # at.spines[ 'right'].set_color(col_grid)
    at.spines['bottom'].set_color(col_frame)
    at.spines[   'top'].set_color(col_frame)
    at.spines[  'left'].set_color(col_frame)
    at.spines[ 'right'].set_color(col_frame)
for at in box_ax:
    at.tick_params(axis = 'both', direction = 'in', color = col_grid, bottom = True, top = True, left = True, right = True, labelsize = fs, length = tl, width = tw)
    # at.spines['bottom'].set_color(col_grid)
    # at.spines[   'top'].set_color(col_grid)
    # at.spines[  'left'].set_color(col_grid)
    # at.spines[ 'right'].set_color(col_grid)
    at.spines['bottom'].set_color(col_frame)
    at.spines[   'top'].set_color(col_frame)
    at.spines[  'left'].set_color(col_frame)
    at.spines[ 'right'].set_color(col_frame)

map_ax[0].set_xlim([Xmax, Xmin])
map_ax[0].set_ylim([Ymin, Ymax])

map_ax[0].set_xlabel(r'$X$ (\si{deg})', fontsize = fs)
map_ax[0].set_ylabel(r'$Y$ (\si{deg})', fontsize = fs)

box_ax[0].set_xlabel(r'Depth ID', fontsize = fs)
box_ax[0].set_ylabel(r'Field ID', fontsize = fs)


# add colorbar
x0, x1 = 1, 0
y0, y1 = 1, 0
for at in map_ax:
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
cax = map_fig.add_axes([x1, y0, 0.05 / nxpanel, y1 - y0])
bar = map_fig.colorbar(map_img, cax = cax)
bar.solids.set_edgecolor('face')
bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-2}})', fontsize = fs)
cax.tick_params(labelsize = fs)

x0, x1 = 1, 0
y0, y1 = 1, 0
for at in box_ax:
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
cax = box_fig.add_axes([x1, y0, 0.05 / nxpanel, y1 - y0])
bar = box_fig.colorbar(box_img, cax = cax)
bar.solids.set_edgecolor('face')
bar.set_label(r'$M$ (\si{M_\odot})', fontsize = fs)
cax.tick_params(labelsize = fs)



# save figures
map_figname = 'fig/' + filename + '_map'
box_figname = 'fig/' + filename + '_box'

map_fig.savefig(map_figname + '.png', format = 'png', dpi =  96, bbox_inches = 'tight')
box_fig.savefig(box_figname + '.png', format = 'png', dpi =  96, bbox_inches = 'tight')
if outputPDF:
    map_fig.savefig(map_figname + '.pdf', format = 'pdf', dpi = 300, bbox_inches = 'tight')
    box_fig.savefig(box_figname + '.pdf', format = 'pdf', dpi = 300, bbox_inches = 'tight')
plt.close('all')




    # for ii in range(nxpanel):
    #     # add colorbar
    #     x0, x1 = 1, 0
    #     y0, y1 = 1, 0
    #     # for at in ax[ii * nypanel:ii * nypanel + (nypanel - 1)]:
    #     for at in ax[ii * nypanel:(ii + 1) * nypanel]:
    #         xl, xr = at.get_position().x0, at.get_position().x1
    #         yb, yt = at.get_position().y0, at.get_position().y1

    #         if x0 > xl:
    #             x0 = xl
    #         if x1 < xr:
    #             x1 = xr
    #         if y0 > yb:
    #             y0 = yb
    #         if y1 < yt:
    #             y1 = yt

    #     cax = fig.add_axes([x0, y1, x1 - x0, 0.05 / nypanel])
    #     if ii == 0:
    #         bar = fig.colorbar(img_Sigma, cax = cax, orientation = 'horizontal')
    #     if ii == 1:
    #         bar = fig.colorbar(img_depth, cax = cax, orientation = 'horizontal')
    #     if ii == 2:
    #         bar = fig.colorbar(img_phase, cax = cax, orientation = 'horizontal')
    #     bar.solids.set_edgecolor('face')
    #     if useDegree:
    #         if ii == 0:
    #             bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-2}})', fontsize = fs)
    #         if ii == 1:
    #             bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-1}.kpc^{-1}})', fontsize = fs)
    #         if ii == 2:
    #             bar.set_label(r'$\Sigma$ (\si{M_\odot.deg^{-1}.km^{-1}.s})', fontsize = fs)
    #     else:
    #         if ii == 0:
    #             bar.set_label(r'$\Sigma$ (\si{M_\odot.kpc^{-2}})', fontsize = fs)
    #         if ii == 1:
    #             bar.set_label(r'$\Sigma$ (\si{M_\odot.kpc^{-2}})', fontsize = fs)
    #         if ii == 2:
    #             bar.set_label(r'$\Sigma$ (\si{M_\odot.kpc^{-1}.km^{-1}.s})', fontsize = fs)
    #     cax.tick_params(labelsize = fs)
    #     cax.xaxis.tick_top()
    #     cax.xaxis.set_label_position('top')


