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


smap_min, smap_max =
vmap_min, vmap_max =


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


init = 50
# last = 100
last = 66


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
    with open('doc/' + filename + '-continue' + '.summary.txt', 'r') as txtfile:
        unit = int(txtfile.readline())
        line = txtfile.readline()
        item = line.split('\t')
        kind = int(item[0])
else:
    kind = None
kind = comm.bcast(kind, root = 0)


cmap = 'plasma'
# cmap = 'cividis'
# cmap = 'inferno'
# cmap = 'viridis'


for fileid in range(init + rank, last + 1, size):
    snapshot = '{:03d}'.format(fileid)
    input_file = 'dat/' + filename + '.subaru' + snapshot + '.h5'

    # read snapshot
    with h5py.File(input_file, 'r') as h5file:
        # read attributes
        Nfield = h5file['/subaru/'].attrs['Nfield'][0]
        Nmodel = h5file['/subaru/'].attrs['Nmodel'][0]
        Nnoise = h5file['/subaru/'].attrs['Nnoise'][0]

        field = h5file['/subaru/'].attrs['field']
        model = h5file['/subaru/'].attrs['model']

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
                xx[kk][ff] = h5file[folder + field[ff] + '-xi']
                yy[kk][ff] = h5file[folder + field[ff] + '-eta']
                vv[kk][ff] = h5file[folder + field[ff] + '-vlos']
                for mm in range(Nmodel):
                    xy_map[kk][ff][mm] = h5file[folder + field[ff] + '-' + model[mm] + '-map']
                    yv_map[kk][ff][mm] = h5file[folder + field[ff] + '-' + model[mm] + '-vel']
                for mm in range(Nnoise):
                    xy_noise[kk][ff][mm] = h5file[folder + field[ff] + '-SN' + '{:02d}'.format(mm + 1) + '-map']
                    yv_noise[kk][ff][mm] = h5file[folder + field[ff] + '-SN' + '{:02d}'.format(mm + 1) + '-vel']
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
                    smap[ii][0] = xy_map[kk][ff][ii] + xy_map[kk][ff][Nmodel + noise]
                    vmap[ii][0] = yv_map[kk][ff][ii] + yv_map[kk][ff][Nmodel + noise]

                    if !specify:
                        smap_min =
                        smap_max =
                        vmap_min =
                        vmap_max = 

                    smap = np.maximum(smap, smap_min)
                    smap = np.minimum(smap, smap_max)
                    vmap = np.maximum(vmap, vmap_min)
                    vmap = np.minimum(vmap, vmap_max)

                    fig = utils.set_figure(nxpanel, nypanel)
                    ax = [0] * nxpanel * nypanel
                    utils.locate_panels(fig, ax, nxpanel, nypanel, True, True)

                    xmin, xmax = xx[0][0], xx[0][nx]
                    ymin, ymax = yy[0][0], yy[0][ny]
                    vmin, vmax = vv[0][0], vv[0][nv]

                    # adjust marker size and etcetra
                    fs = fs_base# / np.sqrt(nypanel)
                    tl = tl_base# / np.sqrt(nypanel)
                    tw = tw_base# / np.sqrt(nypanel)

                    for jj in range(nypanel):
                        img_Sigma = ax[              jj].imshow(xy_map[jj].T, extent = [xmin, xmax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = Sigma_min, vmax = Sigma_max), cmap = cmap_Sigma, aspect = 'auto', rasterized = True)
                        img_phase = ax[2 * nypanel + jj].imshow(yv_map[jj]  , extent = [vmin, vmax, ymin, ymax], origin = 'lower', interpolation = 'none', norm = LogNorm(vmin = phase_min, vmax = phase_max), cmap = cmap_phase, aspect = 'auto', rasterized = True)
