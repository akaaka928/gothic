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
    argparser.add_argument('-c', '--compare',
                           action='store_true',
                           help='Whether to compare with V100 SXM2 (Pascal mode)')
    argparser.add_argument('-p', '--pdf',
                           action='store_true',
                           help='Whether to output figures in PDF format')
    return argparser.parse_args()

# catch optional input arguments
args = get_option()
compare = args.compare
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
model = '1410/log/' + particle + '_'
N001k = pandas.read_csv(model + tag[ 0] + series + str(num[ 0]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N002k = pandas.read_csv(model + tag[ 1] + series + str(num[ 1]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N004k = pandas.read_csv(model + tag[ 2] + series + str(num[ 2]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N008k = pandas.read_csv(model + tag[ 3] + series + str(num[ 3]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N016k = pandas.read_csv(model + tag[ 4] + series + str(num[ 4]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N032k = pandas.read_csv(model + tag[ 5] + series + str(num[ 5]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N064k = pandas.read_csv(model + tag[ 6] + series + str(num[ 6]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N128k = pandas.read_csv(model + tag[ 7] + series + str(num[ 7]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N256k = pandas.read_csv(model + tag[ 8] + series + str(num[ 8]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N512k = pandas.read_csv(model + tag[ 9] + series + str(num[ 9]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N001M = pandas.read_csv(model + tag[10] + series + str(num[10]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N002M = pandas.read_csv(model + tag[11] + series + str(num[11]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N004M = pandas.read_csv(model + tag[12] + series + str(num[12]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N008M = pandas.read_csv(model + tag[13] + series + str(num[13]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
N016M = pandas.read_csv(model + tag[14] + series + str(num[14]) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')


# pick up representative values
best = [min(N001k['elapsed_step']), min(N002k['elapsed_step']), min(N004k['elapsed_step']), min(N008k['elapsed_step']), min(N016k['elapsed_step']), min(N032k['elapsed_step']), min(N064k['elapsed_step']), min(N128k['elapsed_step']), min(N256k['elapsed_step']), min(N512k['elapsed_step']), min(N001M['elapsed_step']), min(N002M['elapsed_step']), min(N004M['elapsed_step']), min(N008M['elapsed_step']), min(N016M['elapsed_step'])]
worst = [max(N001k['elapsed_step']), max(N002k['elapsed_step']), max(N004k['elapsed_step']), max(N008k['elapsed_step']), max(N016k['elapsed_step']), max(N032k['elapsed_step']), max(N064k['elapsed_step']), max(N128k['elapsed_step']), max(N256k['elapsed_step']), max(N512k['elapsed_step']), max(N001M['elapsed_step']), max(N002M['elapsed_step']), max(N004M['elapsed_step']), max(N008M['elapsed_step']), max(N016M['elapsed_step'])]
median = [numpy.median(N001k['elapsed_step']), numpy.median(N002k['elapsed_step']), numpy.median(N004k['elapsed_step']), numpy.median(N008k['elapsed_step']), numpy.median(N016k['elapsed_step']), numpy.median(N032k['elapsed_step']), numpy.median(N064k['elapsed_step']), numpy.median(N128k['elapsed_step']), numpy.median(N256k['elapsed_step']), numpy.median(N512k['elapsed_step']), numpy.median(N001M['elapsed_step']), numpy.median(N002M['elapsed_step']), numpy.median(N004M['elapsed_step']), numpy.median(N008M['elapsed_step']), numpy.median(N016M['elapsed_step'])]
p1sigma = [numpy.percentile(N001k['elapsed_step'], 84.13), numpy.percentile(N002k['elapsed_step'], 84.13), numpy.percentile(N004k['elapsed_step'], 84.13), numpy.percentile(N008k['elapsed_step'], 84.13), numpy.percentile(N016k['elapsed_step'], 84.13), numpy.percentile(N032k['elapsed_step'], 84.13), numpy.percentile(N064k['elapsed_step'], 84.13), numpy.percentile(N128k['elapsed_step'], 84.13), numpy.percentile(N256k['elapsed_step'], 84.13), numpy.percentile(N512k['elapsed_step'], 84.13), numpy.percentile(N001M['elapsed_step'], 84.13), numpy.percentile(N002M['elapsed_step'], 84.13), numpy.percentile(N004M['elapsed_step'], 84.13), numpy.percentile(N008M['elapsed_step'], 84.13), numpy.percentile(N016M['elapsed_step'], 84.13)]
m1sigma = [numpy.percentile(N001k['elapsed_step'], 15.87), numpy.percentile(N002k['elapsed_step'], 15.87), numpy.percentile(N004k['elapsed_step'], 15.87), numpy.percentile(N008k['elapsed_step'], 15.87), numpy.percentile(N016k['elapsed_step'], 15.87), numpy.percentile(N032k['elapsed_step'], 15.87), numpy.percentile(N064k['elapsed_step'], 15.87), numpy.percentile(N128k['elapsed_step'], 15.87), numpy.percentile(N256k['elapsed_step'], 15.87), numpy.percentile(N512k['elapsed_step'], 15.87), numpy.percentile(N001M['elapsed_step'], 15.87), numpy.percentile(N002M['elapsed_step'], 15.87), numpy.percentile(N004M['elapsed_step'], 15.87), numpy.percentile(N008M['elapsed_step'], 15.87), numpy.percentile(N016M['elapsed_step'], 15.87)]
# upper = [numpy.percentile(N001k['elapsed_step'], 75.0), numpy.percentile(N002k['elapsed_step'], 75.0), numpy.percentile(N004k['elapsed_step'], 75.0), numpy.percentile(N008k['elapsed_step'], 75.0), numpy.percentile(N016k['elapsed_step'], 75.0), numpy.percentile(N032k['elapsed_step'], 75.0), numpy.percentile(N064k['elapsed_step'], 75.0), numpy.percentile(N128k['elapsed_step'], 75.0), numpy.percentile(N256k['elapsed_step'], 75.0), numpy.percentile(N512k['elapsed_step'], 75.0), numpy.percentile(N001M['elapsed_step'], 75.0), numpy.percentile(N002M['elapsed_step'], 75.0), numpy.percentile(N004M['elapsed_step'], 75.0), numpy.percentile(N008M['elapsed_step'], 75.0), numpy.percentile(N016M['elapsed_step'], 75.0)]
# lower = [numpy.percentile(N001k['elapsed_step'], 25.0), numpy.percentile(N002k['elapsed_step'], 25.0), numpy.percentile(N004k['elapsed_step'], 25.0), numpy.percentile(N008k['elapsed_step'], 25.0), numpy.percentile(N016k['elapsed_step'], 25.0), numpy.percentile(N032k['elapsed_step'], 25.0), numpy.percentile(N064k['elapsed_step'], 25.0), numpy.percentile(N128k['elapsed_step'], 25.0), numpy.percentile(N256k['elapsed_step'], 25.0), numpy.percentile(N512k['elapsed_step'], 25.0), numpy.percentile(N001M['elapsed_step'], 25.0), numpy.percentile(N002M['elapsed_step'], 25.0), numpy.percentile(N004M['elapsed_step'], 25.0), numpy.percentile(N008M['elapsed_step'], 25.0), numpy.percentile(N016M['elapsed_step'], 25.0)]


fig_box, ax_box, fs_box, ms_box, lw_box = util.set_figure()
if compare:
    v100 = pandas.read_csv('v100_sxm2.csv')
    ax_box[0].plot(v100['Ntot'], v100['elapsed_step'], linestyle = ls[1], linewidth = lw_box, color = col[0], label = 'V100 (SXM2)')

    fig_ratio, ax_ratio, fs_ratio, ms_ratio, lw_ratio = util.set_figure()
    ax_ratio[0].axhline(1555.0 / 900.0, linestyle = ls[2], linewidth = lw_ratio, color = col[0], label = 'memory bandwidth')
    ax_ratio[0].axhline(((64 * 108) * 1.410 * 2) / ((64 * 80) * 1.530 * 2), linestyle = ls[1], linewidth = lw_ratio, color = col[0], label = 'peak performance')
    ax_ratio[0].plot(num, v100['elapsed_step'][v100['Ntot'] <= 16777216] / best, linestyle = ls[0], linewidth = lw_ratio, color = col[0], label = 'measured')
    ax_ratio[0].set_xlabel(r'$N$', fontsize = fs_ratio)
    ax_ratio[0].set_ylabel('Speed up from V100 (SXM2)', fontsize = fs_ratio)
    ax_ratio[0].set_xlim(util.scale_axis(min(num), max(num), True))
    ax_ratio[0].semilogx()
    ax_ratio[0].grid()
    handles, labels = ax_ratio[0].get_legend_handles_labels()
    ax_ratio[0].legend(handles[::-1], labels[::-1], numpoints = 1, loc = 'best', fontsize = fs_ratio)
    fig_ratio.savefig('fig/' + particle + filename + '_speedup' + '.png', format = 'png', dpi = 100, bbox_inches = 'tight')
    if outputPDF:
        fig_ratio.savefig('fig/' + particle + filename + '_speedup' + '.pdf', format = 'pdf', bbox_inches = 'tight')

ax_box[0].plot(num, median, linestyle = ls[2], linewidth = lw_box, color = col[2], label = 'median')
ax_box[0].plot(num, best, linestyle = ls[0], linewidth = lw_box, color = col[1], label = 'fastest')
style_box = dict(lw = 0.75 * lw_box)
style_whi = dict(markersize = ms_box)
style_med = dict(color = col[2], lw = 0.5 * lw_box)
style_fli = dict(marker = 'x', markersize = 0.5 * ms_box)
ax_box[0].boxplot([N001k['elapsed_step'], N002k['elapsed_step'], N004k['elapsed_step'], N008k['elapsed_step'], N016k['elapsed_step'], N032k['elapsed_step'], N064k['elapsed_step'], N128k['elapsed_step'], N256k['elapsed_step'], N512k['elapsed_step'], N001M['elapsed_step'], N002M['elapsed_step'], N004M['elapsed_step'], N008M['elapsed_step'], N016M['elapsed_step']], positions = num, widths = numpy.array(num) * 0.25, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)

ax_box[0].set_xlabel(r'$N$', fontsize = fs_box)
ax_box[0].set_ylabel(r'$t_\mathrm{step}$~(\si{s})', fontsize = fs_box)
ax_box[0].set_xlim(util.scale_axis(min(num), max(num), True))

ax_box[0].loglog()
ax_box[0].grid()

handles, labels = ax_box[0].get_legend_handles_labels()
ax_box[0].legend(handles[::-1], labels[::-1], numpoints = 1, loc = 'best', fontsize = fs_box)



# save figures
fig_box.savefig('fig/' + particle + filename + '_box' + '.png', format = 'png', dpi = 100, bbox_inches = 'tight')
if outputPDF:
    fig_box.savefig('fig/' + particle + filename + '_box' + '.pdf', format = 'pdf', bbox_inches = 'tight')

pyplot.close('all')
