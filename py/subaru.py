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
    argparser.add_argument('-p', '--pdf',
                           action='store_true',
                           help='Whether to output figures in PDF format')
    argparser.add_argument('-m', '--monochrome',
                           action='store_true',
                           help='Whether to draw in monochrome')
    return argparser.parse_args()

# catch optional input arguments
args = get_option()
filename = args.files
outputPDF = args.pdf
monochrome = args.monochrome


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

        nxpanel = Nmodel
        nypanel = 2 # w/o or w/z noise
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
