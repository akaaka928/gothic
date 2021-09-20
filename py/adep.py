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
    argparser.add_argument('-a', '--allvers',
                           action='store_true',
                           help='Whether to compare with all earlier versions (from Fermi to Volta)')
    argparser.add_argument('-f', '--frequency',
                           action='store_true',
                           help='Whether to compare with multiple frequencies')
    argparser.add_argument('-p', '--pdf',
                           action='store_true',
                           help='Whether to output figures in PDF format')
    return argparser.parse_args()

# catch optional input arguments
args = get_option()
compare = args.compare
compare_all = args.allvers
compare_freq = args.frequency
outputPDF = args.pdf
if compare_freq:
    compare = False
if not compare:
    compare_all = False

filename = '_adep' if not compare_all else '_adep_allgen'
if compare_freq:
    filename = '_adep_fdep'

# embed fonts
pyplot.rcParams['ps.useafm'] = True
pyplot.rcParams['pdf.use14corefonts'] = True
pyplot.rcParams['text.usetex'] = True

# use packages
pyplot.rcParams['text.latex.preamble'] = r'\usepackage{physics,siunitx}'

# specify direction of ticks
pyplot.rcParams['xtick.direction'] = 'in'
pyplot.rcParams['ytick.direction'] = 'in'


acc = [5.0000000000000e-01, 2.5000000000000e-01, 1.2500000000000e-01, 6.2500000000000e-02, 3.1250000000000e-02, 1.5625000000000e-02, 7.8125000000000e-03, 3.9062500000000e-03, 1.9531250000000e-03, 9.7656250000000e-04, 4.8828125000000e-04, 2.4414062500000e-04, 1.2207031250000e-04, 6.1035156250000e-05, 3.0517578125000e-05, 1.5258789062500e-05, 7.6293945312500e-06, 3.8146972656250e-06, 1.9073486328125e-06, 9.5367431640625e-07]


# read measured data
particle = 'm31'
series = '.acceleration.block.cc80.'
log = '.time.log'
# measured results with --gpu-freq=1410
boost_freq=1410.0
f1410 = pandas.read_csv('1410/log/' + particle + series + str(8388608) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
f1410_a00 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 0]]
f1410_a01 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 1]]
f1410_a02 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 2]]
f1410_a03 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 3]]
f1410_a04 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 4]]
f1410_a05 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 5]]
f1410_a06 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 6]]
f1410_a07 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 7]]
f1410_a08 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 8]]
f1410_a09 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[ 9]]
f1410_a10 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[10]]
f1410_a11 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[11]]
f1410_a12 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[12]]
f1410_a13 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[13]]
f1410_a14 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[14]]
f1410_a15 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[15]]
f1410_a16 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[16]]
f1410_a17 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[17]]
f1410_a18 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[18]]
f1410_a19 = f1410['elapsed_step'][f1410['accuracy_parameter'] == acc[19]]

if compare_freq:
    # measured results with --gpu-freq=705
    half_freq=705.0
    f0705 = pandas.read_csv('0705/log/' + particle + series + str(8388608) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
    f0705_a00 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 0]]
    f0705_a01 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 1]]
    f0705_a02 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 2]]
    f0705_a03 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 3]]
    f0705_a04 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 4]]
    f0705_a05 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 5]]
    f0705_a06 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 6]]
    f0705_a07 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 7]]
    f0705_a08 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 8]]
    f0705_a09 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[ 9]]
    f0705_a10 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[10]]
    f0705_a11 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[11]]
    f0705_a12 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[12]]
    f0705_a13 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[13]]
    f0705_a14 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[14]]
    f0705_a15 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[15]]
    f0705_a16 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[16]]
    f0705_a17 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[17]]
    f0705_a18 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[18]]
    f0705_a19 = f0705['elapsed_step'][f0705['accuracy_parameter'] == acc[19]]

    f0705_best = numpy.array([min(f0705_a00), min(f0705_a01), min(f0705_a02), min(f0705_a03), min(f0705_a04), min(f0705_a05), min(f0705_a06), min(f0705_a07), min(f0705_a08), min(f0705_a09), min(f0705_a10), min(f0705_a11), min(f0705_a12), min(f0705_a13), min(f0705_a14), min(f0705_a15), min(f0705_a16), min(f0705_a17), min(f0705_a18), min(f0705_a19)])
    f0705_worst = numpy.array([max(f0705_a00), max(f0705_a01), max(f0705_a02), max(f0705_a03), max(f0705_a04), max(f0705_a05), max(f0705_a06), max(f0705_a07), max(f0705_a08), max(f0705_a09), max(f0705_a10), max(f0705_a11), max(f0705_a12), max(f0705_a13), max(f0705_a14), max(f0705_a15), max(f0705_a16), max(f0705_a17), max(f0705_a18), max(f0705_a19)])
    f0705_median = numpy.array([numpy.median(f0705_a00), numpy.median(f0705_a01), numpy.median(f0705_a02), numpy.median(f0705_a03), numpy.median(f0705_a04), numpy.median(f0705_a05), numpy.median(f0705_a06), numpy.median(f0705_a07), numpy.median(f0705_a08), numpy.median(f0705_a09), numpy.median(f0705_a10), numpy.median(f0705_a11), numpy.median(f0705_a12), numpy.median(f0705_a13), numpy.median(f0705_a14), numpy.median(f0705_a15), numpy.median(f0705_a16), numpy.median(f0705_a17), numpy.median(f0705_a18), numpy.median(f0705_a19)])



# pick up representative values
best = numpy.array([min(f1410_a00), min(f1410_a01), min(f1410_a02), min(f1410_a03), min(f1410_a04), min(f1410_a05), min(f1410_a06), min(f1410_a07), min(f1410_a08), min(f1410_a09), min(f1410_a10), min(f1410_a11), min(f1410_a12), min(f1410_a13), min(f1410_a14), min(f1410_a15), min(f1410_a16), min(f1410_a17), min(f1410_a18), min(f1410_a19)])
worst = numpy.array([max(f1410_a00), max(f1410_a01), max(f1410_a02), max(f1410_a03), max(f1410_a04), max(f1410_a05), max(f1410_a06), max(f1410_a07), max(f1410_a08), max(f1410_a09), max(f1410_a10), max(f1410_a11), max(f1410_a12), max(f1410_a13), max(f1410_a14), max(f1410_a15), max(f1410_a16), max(f1410_a17), max(f1410_a18), max(f1410_a19)])
median = numpy.array([numpy.median(f1410_a00), numpy.median(f1410_a01), numpy.median(f1410_a02), numpy.median(f1410_a03), numpy.median(f1410_a04), numpy.median(f1410_a05), numpy.median(f1410_a06), numpy.median(f1410_a07), numpy.median(f1410_a08), numpy.median(f1410_a09), numpy.median(f1410_a10), numpy.median(f1410_a11), numpy.median(f1410_a12), numpy.median(f1410_a13), numpy.median(f1410_a14), numpy.median(f1410_a15), numpy.median(f1410_a16), numpy.median(f1410_a17), numpy.median(f1410_a18), numpy.median(f1410_a19)])
p1sigma = numpy.array([numpy.percentile(f1410_a00, 84.13), numpy.percentile(f1410_a01, 84.13), numpy.percentile(f1410_a02, 84.13), numpy.percentile(f1410_a03, 84.13), numpy.percentile(f1410_a04, 84.13), numpy.percentile(f1410_a05, 84.13), numpy.percentile(f1410_a06, 84.13), numpy.percentile(f1410_a07, 84.13), numpy.percentile(f1410_a08, 84.13), numpy.percentile(f1410_a09, 84.13), numpy.percentile(f1410_a10, 84.13), numpy.percentile(f1410_a11, 84.13), numpy.percentile(f1410_a12, 84.13), numpy.percentile(f1410_a13, 84.13), numpy.percentile(f1410_a14, 84.13), numpy.percentile(f1410_a15, 84.13), numpy.percentile(f1410_a16, 84.13), numpy.percentile(f1410_a17, 84.13), numpy.percentile(f1410_a18, 84.13), numpy.percentile(f1410_a19, 84.13)])
m1sigma = numpy.array([numpy.percentile(f1410_a00, 15.87), numpy.percentile(f1410_a01, 15.87), numpy.percentile(f1410_a02, 15.87), numpy.percentile(f1410_a03, 15.87), numpy.percentile(f1410_a04, 15.87), numpy.percentile(f1410_a05, 15.87), numpy.percentile(f1410_a06, 15.87), numpy.percentile(f1410_a07, 15.87), numpy.percentile(f1410_a08, 15.87), numpy.percentile(f1410_a09, 15.87), numpy.percentile(f1410_a10, 15.87), numpy.percentile(f1410_a11, 15.87), numpy.percentile(f1410_a12, 15.87), numpy.percentile(f1410_a13, 15.87), numpy.percentile(f1410_a14, 15.87), numpy.percentile(f1410_a15, 15.87), numpy.percentile(f1410_a16, 15.87), numpy.percentile(f1410_a17, 15.87), numpy.percentile(f1410_a18, 15.87), numpy.percentile(f1410_a19, 15.87)])


fig_box, ax_box, fs_box, ms_box, lw_box = util.set_figure()
if compare:
    if compare_all:
        m2090 = pandas.read_csv('past/' + particle + '.acceleration.block.' + 'cc20.' + str(8388608) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
        ax_box[0].plot(m2090['accuracy_parameter'], m2090['elapsed_step'], linestyle = ls[6 % Nls], linewidth = lw_box, color = col[7 % Ncol], label = 'Tesla M2090')
        k20x = pandas.read_csv('past/' + particle + '.acceleration.block.' + 'cc35.' + str(8388608) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
        ax_box[0].plot(k20x['accuracy_parameter'], k20x['elapsed_step'], linestyle = ls[5 % Nls], linewidth = lw_box, color = col[3 % Ncol], label = 'Tesla K20X')
        maxwell = pandas.read_csv('past/' + particle + '.acceleration.block.' + 'cc52.' + str(8388608) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
        ax_box[0].plot(maxwell['accuracy_parameter'], maxwell['elapsed_step'], linestyle = ls[4 % Nls], linewidth = lw_box, color = col[5 % Ncol], label = 'GeForce GTX TITAN X')
    p100 = pandas.read_csv('past/' + particle + '.acceleration.block.' + 'cc60.' + str(8388608) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
    ax_box[0].plot(p100['accuracy_parameter'], p100['elapsed_step'], linestyle = ls[3 % Nls], linewidth = lw_box, color = col[4 % Ncol], label = 'Tesla P100 (SXM)')
    v100 = pandas.read_csv('past/' + particle + '.acceleration.block.' + 'cc70_60.' + str(8388608) + log, names = ('accuracy_parameter', 'elapsed_total', 'step', 'elapsed_step'), delimiter = '\t')
    ax_box[0].plot(v100['accuracy_parameter'], v100['elapsed_step'], linestyle = ls[1 % Nls], linewidth = lw_box, color = col[6 % Ncol], label = 'Tesla V100 (SXM2)')

    fig_ratio, ax_ratio, fs_ratio, ms_ratio, lw_ratio = util.set_figure()
    ax_ratio[0].axhline(1555.0 / 900.0, linestyle = ls[2], linewidth = lw_ratio, color = col[0], label = 'memory bandwidth')
    ax_ratio[0].axhline(((64 * 108) * 1.410 * 2) / ((64 * 80) * 1.530 * 2), linestyle = ls[1], linewidth = lw_ratio, color = col[0], label = 'peak performance')
    ax_ratio[0].plot(acc, v100['elapsed_step'] / best, linestyle = ls[0], linewidth = lw_ratio, color = col[0], label = 'measured')
    ax_ratio[0].set_xlabel(r'$\varDelta_\mathrm{acc}$', fontsize = fs_ratio)
    ax_ratio[0].set_ylabel('Speed up from V100 (SXM2)', fontsize = fs_ratio)
    ax_ratio[0].set_xlim(util.scale_axis(min(acc), max(acc), True))
    ax_ratio[0].semilogx()
    ax_ratio[0].grid()
    handles, labels = ax_ratio[0].get_legend_handles_labels()
    ax_ratio[0].legend(handles[::-1], labels[::-1], numpoints = 1, loc = 'best', fontsize = fs_ratio)
    fig_ratio.savefig('fig/' + particle + filename + '_speedup' + '.png', format = 'png', dpi = 100, bbox_inches = 'tight')
    if outputPDF:
        fig_ratio.savefig('fig/' + particle + filename + '_speedup' + '.pdf', format = 'pdf', bbox_inches = 'tight')

if compare_freq:
    style_cap = dict(color = col[2])
    style_box = dict(color = col[2], lw = 0.75 * lw_box)
    style_whi = dict(color = col[2], markersize = ms_box)
    style_med = dict(color = col[2], lw = 0.5 * lw_box)
    style_fli = dict(markeredgecolor = col[2], marker = 'x', markersize = 0.5 * ms_box)
    ax_box[0].boxplot([f0705_a00, f0705_a01, f0705_a02, f0705_a03, f0705_a04, f0705_a05, f0705_a06, f0705_a07, f0705_a08, f0705_a09, f0705_a10, f0705_a11, f0705_a12, f0705_a13, f0705_a14, f0705_a15, f0705_a16, f0705_a17, f0705_a18, f0705_a19], positions = acc, widths = numpy.array(acc) * 0.25, capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)
    style_cap = dict(color = col[1])
    style_box = dict(color = col[1], lw = 0.75 * lw_box)
    style_whi = dict(color = col[1], markersize = ms_box)
    style_med = dict(color = col[1], lw = 0.5 * lw_box)
    style_fli = dict(markeredgecolor = col[1], marker = 'x', markersize = 0.5 * ms_box)
    ax_box[0].boxplot([f1410_a00, f1410_a01, f1410_a02, f1410_a03, f1410_a04, f1410_a05, f1410_a06, f1410_a07, f1410_a08, f1410_a09, f1410_a10, f1410_a11, f1410_a12, f1410_a13, f1410_a14, f1410_a15, f1410_a16, f1410_a17, f1410_a18, f1410_a19], positions = acc, widths = numpy.array(acc) * 0.25, capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)

if not compare_freq:
    ax_box[0].plot(acc, median, linestyle = ls[2], linewidth = lw_box, color = col[2], label = 'median' if not compare else 'A100 (PCIe), median')
    ax_box[0].plot(acc, best, linestyle = ls[0], linewidth = lw_box, color = col[1], label = 'fastest' if not compare else 'A100 (PCIe), fastest')
    style_box = dict(lw = 0.75 * lw_box)
    style_whi = dict(markersize = ms_box)
    style_med = dict(color = col[2], lw = 0.5 * lw_box)
    style_fli = dict(marker = 'x', markersize = 0.5 * ms_box)
    ax_box[0].boxplot([f1410_a00, f1410_a01, f1410_a02, f1410_a03, f1410_a04, f1410_a05, f1410_a06, f1410_a07, f1410_a08, f1410_a09, f1410_a10, f1410_a11, f1410_a12, f1410_a13, f1410_a14, f1410_a15, f1410_a16, f1410_a17, f1410_a18, f1410_a19], positions = acc, widths = numpy.array(acc) * 0.25, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)

ax_box[0].set_xlabel(r'$\varDelta_\mathrm{acc}$', fontsize = fs_box)
ax_box[0].set_ylabel(r'$t_\mathrm{step}$~(\si{s})', fontsize = fs_box)
ax_box[0].set_xlim(util.scale_axis(min(acc), max(acc), True))

ax_box[0].loglog()
ax_box[0].grid()

if not compare_freq:
    handles, labels = ax_box[0].get_legend_handles_labels()
    ax_box[0].legend(handles[::-1], labels[::-1], numpoints = 1, loc = 'best', fontsize = fs_box)



# save figures
fig_box.savefig('fig/' + particle + filename + '_box' + '.png', format = 'png', dpi = 100, bbox_inches = 'tight')
if outputPDF:
    fig_box.savefig('fig/' + particle + filename + '_box' + '.pdf', format = 'pdf', bbox_inches = 'tight')

pyplot.close('all')



# def compare(dat0, dat1, dat2, lab0, lab1, lab2, freq0, freq1, freq2, name):
def compare(dat0, dat1, lab0, lab1, freq0, freq1, name):
    fig, ax, fs, ms, lw = util.set_figure()

    style_cap = dict(color = col[0])
    style_box = dict(color = col[0], lw = lw)
    style_whi = dict(color = col[0], markersize = ms)
    style_med = dict(color = col[0], lw = lw)
    style_fli = dict(markeredgecolor = col[0], marker = 'x', markersize = ms)
    # ax[0].boxplot([dat0, dat1 * freq1 / freq0, dat2 * freq2 / freq0], labels = [lab0, lab1, lab2], capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)
    ax[0].boxplot([dat0, dat1 * freq1 / freq0], labels = [lab0, lab1], capprops = style_cap, boxprops = style_box, whiskerprops = style_whi, medianprops = style_med, flierprops = style_fli)

    if freq0 != 1410.0:
        print('freq0 =', freq0, 'while freq0 = 1410.0 is assumed in this routine')
    ax[0].set_ylabel(r'$\pqty{f_\mathrm{GPU} / 1410} \times t_\mathrm{step}$~(\si{s})', fontsize = fs)
    ax[0].grid()

    # save figures
    fig.savefig('fig/' + particle + '_box_fdep_acc' + name + '.png', format = 'png', dpi = 100, bbox_inches = 'tight')
    if outputPDF:
        fig.savefig('fig/' + particle + '_box_fdep_acc' + name + '.pdf', format = 'pdf', bbox_inches = 'tight')

    pyplot.close('all')


if compare_freq:
    tag = ['2^-01', '2^-02', '2^-03', '2^-04', '2^-05', '2^-06', '2^-07', '2^-08', '2^-09', '2^-10', '2^-11', '2^-12', '2^-13', '2^-14', '2^-15', '2^-16', '2^-17', '2^-18', '2^-19', '2^-20']
    list_boost = [f1410_a00, f1410_a01, f1410_a02, f1410_a03, f1410_a04, f1410_a05, f1410_a06, f1410_a07, f1410_a08, f1410_a09, f1410_a10, f1410_a11, f1410_a12, f1410_a13, f1410_a14, f1410_a15, f1410_a16, f1410_a17, f1410_a18, f1410_a19]
    list_half = [f0705_a00, f0705_a01, f0705_a02, f0705_a03, f0705_a04, f0705_a05, f0705_a06, f0705_a07, f0705_a08, f0705_a09, f0705_a10, f0705_a11, f0705_a12, f0705_a13, f0705_a14, f0705_a15, f0705_a16, f0705_a17, f0705_a18, f0705_a19]

    for ii in range(len(acc)):
        # compare(list_boost[ii], list_half[ii], list_default[ii], r'$f_\mathrm{GPU} = \SI{1410}{MHz}$', r'$f_\mathrm{GPU} = \SI{705}{MHz}$', 'default', boost_freq, half_freq, default_freq, tag[ii])
        compare(list_boost[ii], list_half[ii], r'$f_\mathrm{GPU} = \SI{1410}{MHz}$', r'$f_\mathrm{GPU} = \SI{705}{MHz}$', boost_freq, half_freq, tag[ii])


    fig_ratio, ax_ratio, fs_ratio, ms_ratio, lw_ratio = util.set_figure()
    ax_ratio[0].axhline(boost_freq / half_freq, linestyle = ls[1], linewidth = lw_ratio, color = col[0], label = 'frequency ratio')
    ax_ratio[0].plot(acc, f0705_worst / worst, linestyle = ls[3], linewidth = lw_ratio, color = col[3], label = 'slowest')
    ax_ratio[0].plot(acc, f0705_median / median, linestyle = ls[2], linewidth = lw_ratio, color = col[2], label = 'median')
    ax_ratio[0].plot(acc, f0705_best / best, linestyle = ls[0], linewidth = lw_ratio, color = col[1], label = 'fastest')
    ax_ratio[0].set_xlabel(r'$\varDelta_\mathrm{acc}$', fontsize = fs_ratio)
    ax_ratio[0].set_ylabel(r'$t_{\mathrm{step}, \SI{705}{MHz}} / t_{\mathrm{step}, \SI{1410}{MHz}}$', fontsize = fs_ratio)
    ax_ratio[0].set_xlim(util.scale_axis(min(acc), max(acc), True))
    ax_ratio[0].semilogx()
    ax_ratio[0].grid()
    handles, labels = ax_ratio[0].get_legend_handles_labels()
    ax_ratio[0].legend(handles[::-1], labels[::-1], numpoints = 1, loc = 'best', fontsize = fs_ratio)
    fig_ratio.savefig('fig/' + particle + filename + '_speedup_freq' + '.png', format = 'png', dpi = 100, bbox_inches = 'tight')
    if outputPDF:
        fig_ratio.savefig('fig/' + particle + filename + '_speedup_freq' + '.pdf', format = 'pdf', bbox_inches = 'tight')
    pyplot.close('all')
