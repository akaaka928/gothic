import numpy
import pandas

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker

import util


Npt, pt = util.set_point_type()
Nls, ls = util.set_line_style()
Ncol, col = util.set_color_palette_for_color_universal_design()


from argparse import ArgumentParser
def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-p', '--pdf',
                           action='store_true',
                           help='Whether to output figures in PDF format')
    return argparser.parse_args()

# catch optional input arguments
args = get_option()
outputPDF = args.pdf
filename = '_ndep'

# embed fonts
pyplot.rcParams['ps.useafm'] = True
pyplot.rcParams['pdf.use14corefonts'] = True
pyplot.rcParams['text.usetex'] = True

# use packages
pyplot.rcParams['text.latex.preamble'] = r'\usepackage{physics,siunitx}'

# specify direction of ticks
pyplot.rcParams['xtick.direction'] = 'in'
pyplot.rcParams['ytick.direction'] = 'in'


# read measured data
particle = 'm31'
num = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216]
tag = ['001k', '002k', '004k', '008k', '016k', '032k', '064k', '128k', '256k', '512k', '001M', '002M', '004M', '008M', '016M']
series = '.acceleration.block.cc80.'
log = '.time.log'
# measured results with --gpu-freq=1410
boost_freq=1410.0
boost = '1410/log/' + particle + '_'
boost_N001k = pandas.read_csv(boost + tag[ 0] + series + str(num[ 0]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N002k = pandas.read_csv(boost + tag[ 1] + series + str(num[ 1]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N004k = pandas.read_csv(boost + tag[ 2] + series + str(num[ 2]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N008k = pandas.read_csv(boost + tag[ 3] + series + str(num[ 3]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N016k = pandas.read_csv(boost + tag[ 4] + series + str(num[ 4]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N032k = pandas.read_csv(boost + tag[ 5] + series + str(num[ 5]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N064k = pandas.read_csv(boost + tag[ 6] + series + str(num[ 6]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N128k = pandas.read_csv(boost + tag[ 7] + series + str(num[ 7]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N256k = pandas.read_csv(boost + tag[ 8] + series + str(num[ 8]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N512k = pandas.read_csv(boost + tag[ 9] + series + str(num[ 9]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N001M = pandas.read_csv(boost + tag[10] + series + str(num[10]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N002M = pandas.read_csv(boost + tag[11] + series + str(num[11]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N004M = pandas.read_csv(boost + tag[12] + series + str(num[12]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N008M = pandas.read_csv(boost + tag[13] + series + str(num[13]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
boost_N016M = pandas.read_csv(boost + tag[14] + series + str(num[14]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
# measured results with --gpu-freq=705
half_freq=705.0
half = '0705/log/' + particle + '_'
half_N001k = pandas.read_csv(half + tag[ 0] + series + str(num[ 0]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N002k = pandas.read_csv(half + tag[ 1] + series + str(num[ 1]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N004k = pandas.read_csv(half + tag[ 2] + series + str(num[ 2]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N008k = pandas.read_csv(half + tag[ 3] + series + str(num[ 3]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N016k = pandas.read_csv(half + tag[ 4] + series + str(num[ 4]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N032k = pandas.read_csv(half + tag[ 5] + series + str(num[ 5]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N064k = pandas.read_csv(half + tag[ 6] + series + str(num[ 6]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N128k = pandas.read_csv(half + tag[ 7] + series + str(num[ 7]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N256k = pandas.read_csv(half + tag[ 8] + series + str(num[ 8]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N512k = pandas.read_csv(half + tag[ 9] + series + str(num[ 9]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N001M = pandas.read_csv(half + tag[10] + series + str(num[10]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N002M = pandas.read_csv(half + tag[11] + series + str(num[11]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N004M = pandas.read_csv(half + tag[12] + series + str(num[12]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N008M = pandas.read_csv(half + tag[13] + series + str(num[13]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
half_N016M = pandas.read_csv(half + tag[14] + series + str(num[14]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
# measured results without --gpu-freq
default_freq=1410.0
default = 'default/log/' + particle + '_'
default_N001k = pandas.read_csv(default + tag[ 0] + series + str(num[ 0]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N002k = pandas.read_csv(default + tag[ 1] + series + str(num[ 1]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N004k = pandas.read_csv(default + tag[ 2] + series + str(num[ 2]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N008k = pandas.read_csv(default + tag[ 3] + series + str(num[ 3]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N016k = pandas.read_csv(default + tag[ 4] + series + str(num[ 4]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N032k = pandas.read_csv(default + tag[ 5] + series + str(num[ 5]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N064k = pandas.read_csv(default + tag[ 6] + series + str(num[ 6]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N128k = pandas.read_csv(default + tag[ 7] + series + str(num[ 7]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N256k = pandas.read_csv(default + tag[ 8] + series + str(num[ 8]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N512k = pandas.read_csv(default + tag[ 9] + series + str(num[ 9]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N001M = pandas.read_csv(default + tag[10] + series + str(num[10]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N002M = pandas.read_csv(default + tag[11] + series + str(num[11]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N004M = pandas.read_csv(default + tag[12] + series + str(num[12]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N008M = pandas.read_csv(default + tag[13] + series + str(num[13]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
default_N016M = pandas.read_csv(default + tag[14] + series + str(num[14]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')


fig_box, ax_box, fs_box, ms_box, lw_box = util.set_figure()
# measured results without --gpu-freq
style_cap = dict(color = col[0])
style_box = dict(color = col[0], lw = 0.75 * lw_box)
style_whi = dict(color = col[0], markersize = ms_box)
style_med = dict(color = col[0], lw = 0.5 * lw_box)
style_fli = dict(markeredgecolor = col[0], marker = 'x', markersize = 0.5 * ms_box)
ax_box[0].boxplot([default_N001k['elapsed_step'], default_N002k['elapsed_step'], default_N004k['elapsed_step'], default_N008k['elapsed_step'], default_N016k['elapsed_step'], default_N032k['elapsed_step'], default_N064k['elapsed_step'], default_N128k['elapsed_step'], default_N256k['elapsed_step'], default_N512k['elapsed_step'], default_N001M['elapsed_step'], default_N002M['elapsed_step'], default_N004M['elapsed_step'], default_N008M['elapsed_step'], default_N016M['elapsed_step']], positions = num, widths = numpy.array(num) * 0.25, capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)
# measured results with --gpu-freq=705
style_cap = dict(color = col[2])
style_box = dict(color = col[2], lw = 0.75 * lw_box)
style_whi = dict(color = col[2], markersize = ms_box)
style_med = dict(color = col[2], lw = 0.5 * lw_box)
style_fli = dict(markeredgecolor = col[2], marker = 'x', markersize = 0.5 * ms_box)
ax_box[0].boxplot([half_N001k['elapsed_step'], half_N002k['elapsed_step'], half_N004k['elapsed_step'], half_N008k['elapsed_step'], half_N016k['elapsed_step'], half_N032k['elapsed_step'], half_N064k['elapsed_step'], half_N128k['elapsed_step'], half_N256k['elapsed_step'], half_N512k['elapsed_step'], half_N001M['elapsed_step'], half_N002M['elapsed_step'], half_N004M['elapsed_step'], half_N008M['elapsed_step'], half_N016M['elapsed_step']], positions = num, widths = numpy.array(num) * 0.25, capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)
# measured results with --gpu-freq=1410
style_cap = dict(color = col[1])
style_box = dict(color = col[1], lw = 0.75 * lw_box)
style_whi = dict(color = col[1], markersize = ms_box)
style_med = dict(color = col[1], lw = 0.5 * lw_box)
style_fli = dict(markeredgecolor = col[1], marker = 'x', markersize = 0.5 * ms_box)
ax_box[0].boxplot([boost_N001k['elapsed_step'], boost_N002k['elapsed_step'], boost_N004k['elapsed_step'], boost_N008k['elapsed_step'], boost_N016k['elapsed_step'], boost_N032k['elapsed_step'], boost_N064k['elapsed_step'], boost_N128k['elapsed_step'], boost_N256k['elapsed_step'], boost_N512k['elapsed_step'], boost_N001M['elapsed_step'], boost_N002M['elapsed_step'], boost_N004M['elapsed_step'], boost_N008M['elapsed_step'], boost_N016M['elapsed_step']], positions = num, widths = numpy.array(num) * 0.25, capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)

ax_box[0].set_xlabel(r'$N$', fontsize = fs_box)
ax_box[0].set_ylabel(r'$t_\mathrm{step}$~(\si{s})', fontsize = fs_box)
ax_box[0].set_xlim(util.scale_axis(min(num), max(num), True))

ax_box[0].loglog()
ax_box[0].grid()

# handles, labels = ax_box[0].get_legend_handles_labels()
# ax_box[0].legend(handles[::-1], labels[::-1], numpoints = 1, loc = 'best', fontsize = fs_box)



# save figures
fig_box.savefig('fig/' + particle + filename + '_box_fdep' + '.png', format = 'png', dpi = 100, bbox_inches = 'tight')
if outputPDF:
    fig_box.savefig('fig/' + particle + filename + '_box_fdep' + '.pdf', format = 'pdf', bbox_inches = 'tight')

pyplot.close('all')



def compare(dat0, dat1, dat2, lab0, lab1, lab2, freq0, freq1, freq2, name):
    fig, ax, fs, ms, lw = util.set_figure()

    style_cap = dict(color = col[0])
    style_box = dict(color = col[0], lw = lw)
    style_whi = dict(color = col[0], markersize = ms)
    style_med = dict(color = col[0], lw = lw)
    style_fli = dict(markeredgecolor = col[0], marker = 'x', markersize = ms)
    ax[0].boxplot([dat0, dat1 * freq1 / freq0, dat2 * freq2 / freq0], labels = [lab0, lab1, lab2], capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)

    if freq0 != 1410.0:
        print('freq0 =', freq0, 'while freq0 = 1410.0 is assumed in this routine')
    ax[0].set_ylabel(r'$\pqty{f_\mathrm{GPU} / 1410} \times t_\mathrm{step}$~(\si{s})', fontsize = fs)
    ax[0].grid()

    # save figures
    fig.savefig('fig/' + particle + '_box_fdep_N' + name + '.png', format = 'png', dpi = 100, bbox_inches = 'tight')
    if outputPDF:
        fig.savefig('fig/' + particle + '_box_fdep_N' + name + '.pdf', format = 'pdf', bbox_inches = 'tight')

    pyplot.close('all')



list_boost = [boost_N001k['elapsed_step'], boost_N002k['elapsed_step'], boost_N004k['elapsed_step'], boost_N008k['elapsed_step'], boost_N016k['elapsed_step'], boost_N032k['elapsed_step'], boost_N064k['elapsed_step'], boost_N128k['elapsed_step'], boost_N256k['elapsed_step'], boost_N512k['elapsed_step'], boost_N001M['elapsed_step'], boost_N002M['elapsed_step'], boost_N004M['elapsed_step'], boost_N008M['elapsed_step'], boost_N016M['elapsed_step']]
list_half = [half_N001k['elapsed_step'], half_N002k['elapsed_step'], half_N004k['elapsed_step'], half_N008k['elapsed_step'], half_N016k['elapsed_step'], half_N032k['elapsed_step'], half_N064k['elapsed_step'], half_N128k['elapsed_step'], half_N256k['elapsed_step'], half_N512k['elapsed_step'], half_N001M['elapsed_step'], half_N002M['elapsed_step'], half_N004M['elapsed_step'], half_N008M['elapsed_step'], half_N016M['elapsed_step']]
list_default = [default_N001k['elapsed_step'], default_N002k['elapsed_step'], default_N004k['elapsed_step'], default_N008k['elapsed_step'], default_N016k['elapsed_step'], default_N032k['elapsed_step'], default_N064k['elapsed_step'], default_N128k['elapsed_step'], default_N256k['elapsed_step'], default_N512k['elapsed_step'], default_N001M['elapsed_step'], default_N002M['elapsed_step'], default_N004M['elapsed_step'], default_N008M['elapsed_step'], default_N016M['elapsed_step']]

for ii in range(len(num)):
    compare(list_boost[ii], list_half[ii], list_default[ii], r'$f_\mathrm{GPU} = \SI{1410}{MHz}$', r'$f_\mathrm{GPU} = \SI{705}{MHz}$', 'default', boost_freq, half_freq, default_freq, tag[ii])
