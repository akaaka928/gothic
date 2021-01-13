import numpy

import matplotlib.pyplot as pyplot
def set_figure(nx = 1, ny = 1, share_xaxis = True, share_yaxis = True, chars = 24.0, dpi = 300.0, inch = 8.0, xscale = 1.0, yscale = 1.0):
    # chars: number of characters; for setting fontsize
    # dpi: dots per inch; for setting resolution
    # inch: size of panel in units of inch (A4 is 8.27 inch * 14.32 inch)
    fig = pyplot.figure(figsize = (inch * xscale * nx, inch * yscale * ny), dpi = dpi)

    # set default size in units of point
    fontsize = int(inch * numpy.sqrt(max(nx * xscale, ny * yscale)) * 72 / chars) # 72 pt = 1 inch
    markersize = inch
    linewidth = inch * 0.25
    ticklength = 6.0 * linewidth

    # set axes
    ax = [0] * nx * ny
    xmin, xmax = 0.0, 1.0
    ymin, ymax = 0.0, 1.0
    xbin = (xmax - xmin) / nx
    ybin = (ymax - ymin) / ny
    xmargin, ymargin = 0, 0

    margin = 0.15
    if not share_yaxis:
        xmin = 0.0
        xbin = 1.0 / nx
        xmargin = xbin * margin
    if not share_xaxis:
        ymin = 0.0
        ybin = 1.0 / ny
        ymargin = ybin * margin

    for ii in range(nx):
        xl = xmin + ii * xbin + xmargin

        for jj in range(ny):
            yl = ymin + jj * ybin + ymargin
            kk = ii * ny + jj
            ax[kk] = fig.add_axes((xl, yl, xbin - 2 * xmargin, ybin - 2 * ymargin))

    for at in ax:
        for axis in ['top', 'bottom', 'left', 'right']:
            at.spines[axis].set_linewidth(linewidth)
        at.tick_params(axis = 'both', direction = 'in', bottom = True, top = True, left = True, right = True, labelsize = fontsize, length = ticklength, width = linewidth)
        at.tick_params(axis = 'x', pad = 0.3 * fontsize)
        at.tick_params(axis = 'both', which = 'minor', direction = 'in', bottom = True, top = True, left = True, right = True, length = 0.5 * ticklength, width = 0.5 * linewidth)
        # grid_linewidth = linewidth

        if share_xaxis:
            at.tick_params(labelbottom = False)

        if share_yaxis:
            at.tick_params(labelleft = False)

    if share_xaxis:
        for ii in range(nx):
            ax[ii * ny].tick_params(labelbottom = True)

    if share_yaxis:
        for jj in range(ny):
            ax[jj].tick_params(labelleft = True)


    return fig, ax, fontsize, markersize, linewidth


def scale_axis(minimum, maximum, logPlt):
    blankVal = 0.2
    if logPlt:
        width = numpy.log10(maximum / minimum)
        blank = width * blankVal / 2
        scale = numpy.power(10.0, blank)
        return minimum / scale, maximum * scale
    else:
        width = maximum - minimum
        blank = width * blankVal / 2
        right = maximum + blank
        left = minimum
        if minimum != 0.0:
            left -= blank
        return left, right


import sys
def add_colorbar(fig, ax, img, width, label, fs, lw, vertical = True, multiple = True):
    if multiple:
        x0, x1 = sys.float_info.max, -sys.float_info.max
        y0, y1 = sys.float_info.max, -sys.float_info.max
        for at in ax:
            xl, xr = at.get_position().x0, at.get_position().x1
            yb, yt = at.get_position().y0, at.get_position().y1
            x0 = min(x0, xl)
            x1 = max(x1, xr)
            y0 = min(y0, yb)
            y1 = max(y1, yt)
    else:
        x0, x1 = ax.get_position().x0, ax.get_position().x1
        y0, y1 = ax.get_position().y0, ax.get_position().y1

    if vertical:
        cax = fig.add_axes([x1, y0, width, y1 - y0])
        bar = fig.colorbar(img, cax = cax)
        bar.solids.set_edgecolor('face')
        bar.set_label(label, fontsize = fs)
        cax.tick_params(axis = 'y', labelsize = fs, left = False, right = True, length = 12.0 * lw, width = lw, labelleft = False, labelright = True)
        cax.tick_params(axis = 'y', which = 'minor', left = False, right = True, length = 6.0 * lw, width = 0.5 * lw)
    else:
        cax = fig.add_axes([x0, y1, x1 - x0, width])
        bar = fig.colorbar(img, cax = cax, orientation = 'horizontal')
        bar.solids.set_edgecolor('face')
        bar.set_label(label, fontsize = fs, labelpad = 0.5 * fs)
        cax.tick_params(axis = 'x', labelsize = fs, bottom = False, top = True, length = 12.0 * lw, width = lw, labelbottom = False, labeltop = True)
        cax.tick_params(axis = 'x', which = 'minor', bottom = False, top = True, length = 6.0 * lw, width = 0.5 * lw)
        cax.xaxis.set_label_position('top')



def set_shared_xlabel(ax, xlabel):
    fig = ax[-1].figure
    fig.canvas.draw()

    # get the corner for all plots
    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax:
        at.set_xlabel('') # remove existing xlabels

        xl, xr = at.get_position().x0, at.get_position().x1

        bboxes, _ = at.xaxis.get_ticklabel_extents(fig.canvas.renderer)
        bboxes = bboxes.inverse_transformed(fig.transFigure)
        yb, yt = bboxes.y0, bboxes.y1

        if x0 > xl:
            x0 = xl
        if x1 < xr:
            x1 = xr
        if y0 > yb:
            y0 = yb
        if y1 < yt:
            y1 = yt

    # set position of label
    ax[-1].set_xlabel(xlabel)
    # ax[-1].xaxis.set_label_coords((x0 + x1) / 2, (y0 + y1) / 2, transform = fig.transFigure)
    ax[-1].xaxis.set_label_coords((x0 + x1) / 2, y0, transform = fig.transFigure)


def set_shared_ylabel(ax, ylabel):
    fig = ax[-1].figure
    fig.canvas.draw()

    # get the corner for all plots
    x0, x1 = 1, 0
    y0, y1 = 1, 0
    for at in ax:
        at.set_ylabel('') # remove existing ylabels

        yb, yt = at.get_position().y0, at.get_position().y1

        bboxes, _ = at.yaxis.get_ticklabel_extents(fig.canvas.renderer)
        bboxes = bboxes.inverse_transformed(fig.transFigure)
        xl, xr = bboxes.x0, bboxes.x1

        if x0 > xl:
            x0 = xl
        if x1 < xr:
            x1 = xr
        if y0 > yb:
            y0 = yb
        if y1 < yt:
            y1 = yt

    # set position of label
    ax[-1].set_ylabel(ylabel)
    ax[-1].yaxis.set_label_coords((x0 + x1) / 2, (y0 + y1) / 2, transform = fig.transFigure)


def set_global_title(fig, title, fontsize, offset):
    y1 = 0
    for at in fig.get_axes():
        yt = at.get_position().y1
        if y1 < yt:
            y1 = yt

    fig.suptitle(title, fontsize = fontsize, y = offset + y1, verticalalignment = 'bottom')


def set_point_type():
    Npt = 5
    pt = ['o', 's', '^', 'D', 'x']

    return Npt, pt


def set_line_style():
    # Nls = 4
    # ls = ['-', '--', ':', '-.']
    Nls = 5
    ls = ['solid', (0, (1, 1)), (0, (5, 5)), (0, (5, 1, 1, 1)), (0, (5, 1, 1, 1, 1, 1, 1, 1))]

    return Nls, ls


def set_color_palette_for_color_universal_design():
    Ncol = 10
    col = [0] * Ncol

    # taken from Model Color Pallete for Color Universal Design ver.4 (pages 7 and 2)
    # conversion using https://hogehoge.tk/tool/number.html
    col[0] = '#000000'# black
    col[1] = '#ff4b00'# red
    col[2] = '#005aff'# blue
    col[3] = '#f6aa00'# orange
    col[4] = '#03af7a'# green
    col[5] = '#4dc4ff'# sky blue
    col[6] = '#804000'# brown
    col[7] = '#990099'# purple
    col[8] = '#fff100'# yellow
    col[9] = '#ff8082'# pink

    return Ncol, col


def set_color_palette_for_monochrome_color():
    Ncol = 4
    col = [0] * Ncol

    # taken from Model Color Pallete for Color Universal Design ver.4 (page 2)
    # conversion using https://hogehoge.tk/tool/number.html
    col[0] = '#000000'# black
    col[1] = '#84919e'# dark gray
    col[2] = '#c8c8cb'# light gray
    col[3] = '#ffffff'# white

    return Ncol, col


# for scientific notation in coutour plot
# based on http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib
def scientific(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


# generate color map
from matplotlib.colors import LinearSegmentedColormap
def truncate_cmap(old, fmin = 0.0, fmax = 1.0, num = -1):
    cmap = pyplot.get_cmap(old)
    if num < 0:
        num = cmap.N
    return LinearSegmentedColormap.from_list('trunc_cmap', cmap(numpy.linspace(fmin, fmax, num)))

def generate_cmap(cols):
    vals = range(len(cols))
    vmax = numpy.ceil(numpy.max(vals))
    color_list = []
    for vv, cc in zip(vals, cols):
        color_list.append((vv / vmax, cc))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)

def generate_cmap_comp(cols, loci):
    if len(cols) != len(loci):
        print('colors and position list need the same number of elements.')
        sys.exit(0)
    vals = range(len(cols))
    vmax = numpy.ceil(numpy.max(vals))
    color_list = []
    for pp, cc in zip(loci, cols):
        color_list.append((pp, cc))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)
