"""
Chi2 primer
-----------
A simple chi-square test statistic which tests how close an observation is
to an expectated value (a model) is::

    Chi2 = (Obs - Exp)**2 / Exp

The Chi2 statistic (fit value) used by MATCH is a bit more complicated and
is based on Poisson probability, but the idea is the same. The "null
hypothesis" (H0) is the assertion that Obs == Exp, such that Chi2 == 0,
i.e., that the observations are perfectly described by the model. The
"alternative hypothesis" (H1) is that the observations do not match the
model. H0 and H1 may be accepted/rejected based on the confidence intervals
in the Chi2 distribution, which is (one degree of freedom, k=1)::

    f(Chi2) = 1 / sqrt(2 * pi * Chi2 * e^Chi2)

For a given value of Chi2, the confidence interval measures the probability
that H0/H1 is false/true. A small confidence interval in this case would
indicate that the observations are consistent with the proposed model.
Confidence interval is just the area under f(Chi2) from 0 to Chi2.

In Gaussian statistics, 68% of observations will be within one standard deviation (1sigma) of the mean, 95% will be within 2sigma, etc. These
percentages are also confidence intervals. sigma can be defined for a Chi2
distribution using the same concept: let 1sigma correspond to the Chi2 at
which the confidence interval (integral of f(Chi2)) equals 68%, etc. For the
k=1 Chi2 distribution, the Chi2 values at different multiples of sigma levels
are,

===== ====
sigma Chi2
===== ====
1     1
2     4
3     9
n     n**2
===== ====

So, two SFHs calculated by MATCH are consistent within 1sigma in terms of
how "good" the solutions are if their fit values are within 1 of each
other, 2sigma if the difference in fits is 4, etc.

"""
from matplotlib import colors, cm, pyplot as plt, ticker
import numpy as np
import os
from scipy import interpolate
import sys
try:
    from lxml import etree
except ImportError:
    from xml.etree import ElementTree as etree

from researchtools.uvregions import config as c, utils as u


#name_list = sys.argv[1:]
name_list = c.NAME_LIST


xlim, ylim = (0.0, 1.6), (0.0, 2.8)  # Data limits
Avlim = ylim[1]  # Avf+dAv (or x+y) may not exceed this total extinction
mode = 'linear'  # 'linear' or 'voronoi'
insideonly = True  # Only consider points within xlim and ylim and Avlim

res = 0.05  # Minimum step size in Av and dAv
npts = 300  # Num of interp points (larger is smoother, voronoi mode only)
nsig1 = 10  # imshow color map limit
nsig2 = 5  # Point/contour color map limit
levels = [1, 2, 3, 4, 5, 10, 20, 40, 80, 160, 320]


# Custom colormap
cdict = {'red':   ((0.00, 0.80, 0.80),
                   (0.20, 1.00, 1.00),
                   (0.40, 0.00, 0.00),
                   (0.60, 0.00, 0.00),
                   (0.80, 0.30, 0.30),
                   (1.00, 0.80, 0.80)),
         'green': ((0.00, 0.00, 0.00),
                   (0.20, 0.50, 0.50),
                   (0.40, 0.80, 0.80),
                   (0.60, 0.40, 0.40),
                   (0.80, 0.00, 0.00),
                   (1.00, 0.00, 0.00)),
         'blue':  ((0.00, 0.00, 0.00),
                   (0.20, 0.00, 0.00),
                   (0.40, 0.00, 0.00),
                   (0.60, 0.80, 0.80),
                   (0.80, 1.00, 1.00),
                   (1.00, 0.80, 0.80))}
cmap = colors.LinearSegmentedColormap('custom', cdict)
cmap.set_over((0.9, 0.5, 0.0))


def sigma_labels(x, pos):
    return '{0:.0f}'.format(x) + r'$\sigma$'


for name in name_list:
    # Get data from xmlfile
    table = u.read_table(c.regxml_file(name))
    x_list, y_list, fit_list = table['A_V'], table['dA_V'], table['fit']

    # Filter out points outside Avlim and dAvlim
    if insideonly:
        k = ((xlim[0] <= x_list) & (x_list <= xlim[1]) &
             (ylim[0] <= y_list) & (y_list <= ylim[1]) &
             (x_list + y_list <= Avlim))
        x_list, y_list, fit_list = x_list[k], y_list[k], fit_list[k]

    # Best fit; note that the best fit could occur at more than 1 point
    imin = fit_list.argmin()
    xmin, ymin, fitmin = x_list[imin], y_list[imin], fit_list[imin]

    # Data ranges for plot window
    x1, x2 = xlim[0]-res/2, xlim[1]+res/2
    y1, y2 = ylim[0]-res/2, ylim[1]+res/2

    # Interpolate grid for imshow
    nsigma_list = np.sqrt(fit_list - fit_list.min())
    if mode == 'voronoi':
        xbins, ybins = np.linspace(x1, x2, npts), np.linspace(y1, y2, npts)
        interp = 'nearest'
    elif mode == 'linear':
        nx = int(np.round((xlim[1]-xlim[0])/res+1))
        ny = int(np.round((ylim[1]-ylim[0])/res+1))
        xbins = np.linspace(xlim[0], xlim[1], nx)
        ybins = np.linspace(ylim[0], ylim[1], ny)
        interp = 'linear'
    xgrid, ygrid = np.meshgrid(xbins, ybins)
    nsigmagrid = interpolate.griddata(np.c_[x_list, y_list], nsigma_list,
                                      (xgrid, ygrid), method=interp)


    # Figure layout, everything in inches
    #
    #     General settings
    aspect = (y2-y1) / (x2-x1)
    ppi = 72.0  # Points per inch; convert font size in points to inches
    label_size, tick_size, text_size = 10, 10, 10  # points, not inches
    lineheight = 1.2  # font size multiplier
    #
    #     Lengths along x-axis, from left to right
    fig_dx = 4.0  # total figure width
    xs1 = 0.05
    yl1_dx = lineheight * label_size / ppi  # ax1 y-axis label
    xs2 = 0.05
    ytl1_dx = 0.2  # ax1 y-axis tick labels
    xs3 = 0.05
    ax2_dx = 0.15  # ax2
    ytl2_dx = 0.4  # ax2 y-axis tick labels
    xs4 = 0.05
    tb_dx = lineheight * label_size / ppi  # text box
    xs5 = 0.05
    #
    #     Derived lengths
    ax1_dx = fig_dx - sum([xs1, yl1_dx, xs2, ytl1_dx, xs3,
                           ax2_dx, ytl2_dx, xs4, tb_dx, xs5])  # ax1
    ax3_dx = ax1_dx  # ax3
    #
    #     Lengths along y-axis, from bottom to top
    ys1 = 0.05
    xl1_dy = lineheight * label_size / ppi  # ax1 x-axis label
    ys2 = 0.05
    xtl1_dy = lineheight * tick_size / ppi  # ax1 x-axis tick labels
    ys3 = 0.05
    ax3_dy = 0.15  # ax3
    xtl3_dy = lineheight * tick_size / ppi  # ax3 x-axis tick labels
    ys4 = 0.05
    #
    #     Derived lengths
    ax1_dy = aspect * ax1_dx  # ax1
    ax2_dy = ax1_dy  # ax2
    fig_dy = sum([ys1, xl1_dy, ys2, xtl1_dy, ax1_dy,
                  ys3, ax3_dy, xtl3_dy, ys4])  # total figure height
    #
    #     Derived locations
    ax1_x0 = sum([xs1, yl1_dx, xs2, ytl1_dx])
    ax1_y0 = sum([ys1, xl1_dy, ys2, xtl1_dy])
    ax2_x0 = sum([ax1_x0, ax1_dx, xs3])
    ax2_y0 = ax1_y0
    ax3_x0 = ax1_x0
    ax3_y0 = sum([ax1_y0, ax1_dy, ys3])
    xl1_xm = ax1_x0 + ax1_dx/2  # Midpoint of ax1 x-axis label
    yl1_ym = ax1_y0 + ax1_dy/2  # Midpoint of ax1 y-axis label
    tb_x0, tb_y0 = ax2_x0+ax2_dx+ytl2_dx+xs4, ax1_y0  # text box

    # Plot Av,dAv surface
    fig = plt.figure(figsize=(fig_dx, fig_dy))
    ax1 = fig.add_axes([ax1_x0/fig_dx, ax1_y0/fig_dy,
                        ax1_dx/fig_dx, ax1_dy/fig_dy])
    ax1.axis([x1, x2, y1, y2])
    ax1.imshow(nsigmagrid, extent=(x1, x2, y1, y2), origin='lower',
               cmap=cm.bone_r, interpolation='nearest', aspect='normal',
               vmin=0, vmax=nsig1)

    # Plot data points and contours
    ax1.scatter(x_list, y_list, marker='o', s=2, edgecolors='none',
                c=nsigma_list, cmap=cmap, vmin=0, vmax=nsig2)
    ax1.plot(xmin, ymin, 'rx')  # Best fit
    ax1.contour(xgrid, ygrid, nsigmagrid, levels, cmap=cmap, vmin=0, vmax=nsig2,
                alpha=0.5)

    # Ticks, grid
    ax1.tick_params(labelsize=tick_size)
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax1.grid(color='w', linestyle=':', alpha=0.3, which='minor')
    ax1.set_axisbelow(True)

    # Labels
    ax1.text(xl1_xm/fig_dx, ys1/fig_dy, r'$A_V$', ha='center',
             size=label_size, transform=fig.transFigure)
    ax1.text(xs1/fig_dx, yl1_ym/fig_dy, r'$dA_V$', va='center',
             size=label_size, rotation='vertical', transform=fig.transFigure)

    # imshow colorbar
    ax2 = fig.add_axes([ax2_x0/fig_dx, ax2_y0/fig_dy,
                        ax2_dx/fig_dx, ax2_dy/fig_dy])
    arr = np.c_[np.linspace(0, nsig1, 100)]
    ax2.imshow(arr, extent=(0, 1, 0, nsig1), origin='lower', cmap=cm.bone_r,
               aspect='normal', vmin=0, vmax=nsig1)
    ax2.tick_params(labelsize=tick_size, bottom=False, top=False,
                    labelbottom=False, labelleft=False, labelright=True)
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(sigma_labels))

    # Point/contour colorbar
    ax3 = fig.add_axes([ax3_x0/fig_dx, ax3_y0/fig_dy,
                        ax3_dx/fig_dx, ax3_dy/fig_dy])
    arr = np.c_[np.linspace(0, nsig2, 100)].transpose()
    ax3.imshow(arr, extent=(0, nsig2, 0, 1), origin='lower', cmap=cmap,
               aspect='normal', vmin=0, vmax=nsig2)
    ax3.tick_params(labelsize=tick_size, left=False, right=False,
                    labelbottom=False, labelleft=False, labeltop=True)
    ax3.xaxis.set_major_formatter(ticker.FuncFormatter(sigma_labels))

    # Summary text
    text = ('Best: Av, dAv, fit = {0:.2f}, {1:.2f}, {2:.2f}'
            .format(xmin, ymin, fitmin))
    ax2.text(tb_x0/fig_dx, tb_y0/fig_dy, text, size=label_size, va='bottom',
             rotation='vertical', transform=fig.transFigure)

    # Labels, write
    out_name = os.path.join(c.regsfh_dir(name), 'fits_{0:s}.pdf'.format(name))
    fig.savefig(out_name, dpi=200)
    plt.close(fig)
