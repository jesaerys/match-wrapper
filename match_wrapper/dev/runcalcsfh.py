"""
Functions for executing multiple runs of the calcsfh command line program.

Included is a wrapper for the calcsfh command line program, and some Av,
dAv search algorithms to control where calcsfh is run in Av, dAv space.
These functions are written very generally and are not intended to be
customized. All customization (logging, file organization, calcsfh
parameters, etc.) should should be done elsewhere, e.g., main().

The required command line arguments depend on the mode of operation::

  {name} pattern {Av0} {dAv0} {d_Av} {d_dAv} {d_Av0} {d_dAv0}
  {name} grid {Av1} {dAv1} {Av2} {dAv2} {nAv} {ndAv}
  {name} list {pointfile}

`name` could be an object or region name being analyzed, and pattern, grid,
and list refer to functions that manage the calcsfh runs. The arguments are
described in their respective functions.

"""
import errno
import numpy as np
import os
import subprocess
import sys
import time
try:
    from lxml import etree
except ImportError:
    from xml.etree import ElementTree as etree


###TEST_CODE_BEGIN
### Enable experimental features in test mode using TEST = True. All new or
### modified parts of the code are wrapped in ###TEST_CODE tags. if-else
### statements are used to control which code is run: the new/test code and
### the original/non-test code occurs in the if clause and else clause,
### respectively.
TEST = False
###TEST_CODE_END


def format_table(table_xml, level=1, indent='  '):
    """XML pretty-print formatter for TABLE elements."""
    table_xml.text = '\n' + level*indent
    coldef_xml, data_xml = table_xml

    coldef_xml.text = '\n' + (level+1)*indent
    for c_xml in coldef_xml:
        c_xml.tail = '\n' + (level+1)*indent
    c_xml.tail = '\n' + (level)*indent
    coldef_xml.tail = '\n' + (level)*indent

    data_xml.text = '\n' + (level+1)*indent
    if len(data_xml):
        for r_xml in data_xml:
            r_xml.tail = '\n' + (level+1)*indent
        r_xml.tail = '\n' + (level)*indent
    data_xml.tail = '\n' + (level-1)*indent
    table_xml.tail = '\n'


def safe_mkdir(path):
    """Create a directory only if it does not exist."""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def calcsfh(photfile, fakefile, IMF, dmod, Av, res, logZ,
            BF, badfrac, CMD, tbins,
            paramfile=None, sfhfile=None, cmdfile=None, terminalfile=None,
            zinc=None, dAv=None, dAvy=None, models=None, path=None):
    """
    Python wrapper for calcsfh in MATCH.

    Most of the function arguments are used to create a calcsfh parameter
    file on the fly before running calcsfh. The remaining arguments are
    used to set up the calcsfh command. Refer to the MATCH README file for
    further information about the parameter file and useful settings for
    the arguments.

    When the calcsfh has finished running, the fit value of the best-fit
    solution is returned along with all content in the .sfh and .sfh.cmd
    output files. All output that calcsfh sends to stdout is piped to
    terminalfile.

    This wrapper has only been tested with MATCH 2.4. The following
    features are not currently supported:

    - Multi-CMD fitting, i.e., obtaining SFHs from photometry for more than
      two filters.
    - Specification of a background CMD.
    - Many other calcsfh flags (e.g., -logterrsig, -mbolerrsig, etc.).

    Parameters
    ----------
    photfile : string
        Path to the input photometry file.
    fakefile : string
        Path to the fake star photometry file.
    IMF : {'Kroupa', 'Salpeter', float}
        IMF slope. 'Kroupa' is equivalent to -1.0, 'Salpeter' is equivalent
        to 1.35. (not -1.35).
    dmod : list or tuple
        Minimum and maximum distance modulus values, (min, max).
    Av : list or tuple
        Minimum and maximum foreground V-band extinction values,
        (min, max).
    res : float
        Resolution for dmod and Av (the same value is used for both).
    logZ : list or tuple
        Minimum, maximum, and resolution values for metallicity (solar=0),
        (min, max, res).
    BF : int or float
        Binary fraction.
    badfrac : list or tuple
        Fraction of bad detections at the top and bottom of the CMD,
        (top, bottom).

    CMD : dict
        The CMD dict contains the following keys:

        V : dict
            A filter dictionary for the bluer ("V") filter. Contains the
            keys 'name' for the filter name, and 'min' and 'max' for the
            bright and faint magnitude limits of the CMD. E.g., ::

              filter1 = {'name': 'WFC475W', 'min': 16.0, 'max': 27.0}.

        I : dict
            Same as V, but for the redder ("I") filter.
        Vres : float
            Resolution (bin size) of the V filter magnitude in the CMD.
        V-I : list or tuple
            Minimum, maximum, and resolution (bin size) values for V-I
            color in the CMD, (min, max, res).
        fake_sm : int or float
            Factor used to improve fake star statistics.
        exclude_gates : list, optional
            If specified (or if not None), all exclude gates in the list
            will be applied to the CMD. A single gate is a list of four
            points::

              exgate = [(color1, mag1), (color2, mag2),
                        (color3, mag3), (color4, mag4)]

            Clockwise vs. counterclockwise order does not matter.
        combine_gates : list, optional
            Same as exclude_gates, but for combine gates.

    tbins : list
        List of time bins, one bin per element. Each element must be a list
        containing at least the young and old age limits of the bin.
        Additional values may be specified after the old age limit, e.g., a
        third number to fix the SFR in a particular bin.
    paramfile : string, optional
        Path to the calcsfh parameter file created based on the function
        arguments. Default is photfile + '.par'.
    sfhfile : string, optional
        Path to the .sfh output file. Default is photfile + '.sfh'.
    cmdfile : string, optional
        Path to the .cmd output file. Default is photfile + '.sfh.cmd'
    terminalfile : string, optional
        Path to the calcsfh terminal output file. Default is
        photfile + '.out'.
    zinc : list or tuple, optional
        Constrain metallicity to increase with time with a list of minimum and
        maximum initial and final metallicity vaues, (initial_min, initial_max,
        final_min, final_max). Default is None (no constraint).
    dAv : float, optional
        Maximum differential extinction. Default is None (no differential
        extinction).
    dAvy : float, optional
        Maximum additional differential extinction for young stars. Default
        is None (no additional differential extinction).
    models : {'PADUA_AGB', etc.}, optional
        Model selection flag. If None, the calcsfh internal default Padua
        2006 models are used.
    path : string, optional
        Path to the calcsfh executable, e.g.,
        '/usr/local/MATCH/match2.4/bin/calcsfh'. If None (default), then the
        'calcsfh' command is resolved using the $PATH variable.

    Returns
    -------
    fit : float
        Fit value of best-fit solution from calcsfh.
    sfh : string
        All content from the output .sfh file (large amount of data).
    cmd : string
        All content from the output .cmd file (large amount of data).

    """
    if IMF == 'Kroupa':
        IMF = '-1.0'
    elif IMF == 'Salpeter':
        IMF = '1.35'
    else:
        IMF = '{0:.2f}'.format(IMF)
    if paramfile is None:
        paramfile = '{0:s}.par'.format(photfile)
    if sfhfile is None:
        sfhfile = '{0:s}.sfh'.format(photfile)
    if cmdfile is None:
        cmdfile = '{0:s}.sfh.cmd'.format(photfile)
    if terminalfile is None:
        terminalfile = '{0:s}.out'.format(photfile)

    # Write the parameter file
    with open(paramfile, 'w') as f:
        # IMF, dmod, Av
        f.write('{0:s} {1:.2f} {2:.2f} {3:.2f} {4:.2f} {5:.2f} {3:.2f}\n'
                .format(IMF, dmod[0], dmod[1], res, Av[0], Av[1]))

        # logZ, zinc
        line = '{0:.1f} {1:.1f} {2:.2f}'.format(*logZ)
        if zinc is not None:
            line = '{0:s} {1:.1f} {2:.1f} {3:.1f} {4:.1f}\n'.format(line, *zinc)
        else:
            line = '{0:s}\n'.format(line)
        f.write(line)

        # binary fraction, bad fractions
        f.write('{0:.2f} {1:.6f} {2:.6f}\n'.format(BF, badfrac[0], badfrac[1]))

        # CMD info
        N = 1
        f.write('{0:d}\n'.format(N))
        f.write('{0:.2f} {1:.2f} {2:.0f} {3:.2f} {4:.2f} {5:s},{6:s}\n'
                .format(CMD['Vres'], CMD['V-I'][2], CMD['fake_sm'],
                        CMD['V-I'][0], CMD['V-I'][1],
                        CMD['V']['name'], CMD['I']['name']))

        # Filter info
        f.write('{0:.2f} {1:.2f} {2:s}\n'
                .format(CMD['V']['min'], CMD['V']['max'], CMD['V']['name']))
        f.write('{0:.2f} {1:.2f} {2:s}\n'
                .format(CMD['I']['min'], CMD['I']['max'], CMD['I']['name']))

        # Gate info
        if CMD.has_key('exclude_gates') and CMD['exclude_gates'] is not None:
            N = len(CMD['exclude_gates'])
            line = '{0:d}'.format(N)
            for gate in CMD['exclude_gates']:
                for x, y in gate:
                    line = ('{0:s} {1:.2f} {2:.2f}'.format(line, x, y))
        else:
            line = '0'
        if CMD.has_key('combine_gates') and CMD['combine_gates'] is not None:
            N = len(CMD['combine_gates'])
            line = '{0:s} {1:d}'.format(line, N)
            for gate in CMD['combine_gates']:
                for x, y in gate:
                    line = ('{0:s} {1:.2f} {2:.2f}'.format(line, x, y))
        else:
            line = '{0:s} 0'.format(line)
        f.write('{0:s}\n'.format(line))

        # Time bins
        N = len(tbins)
        f.write('{0:d}\n'.format(N))
        for tbin in tbins:
            if len(tbin) == 3:
                line = '{0:.2f} {1:.2f} {2:.4e}\n'.format(*tbin)
            else:
                line = '{0:.2f} {1:.2f}\n'.format(*tbin)
            f.write(line)


    # Construct the calcsfh command
    if path is None:
        path = 'calcsfh'
    command = ('{0:s} {1:s} {2:s} {3:s} {4:s}'
               .format(path, paramfile, photfile, fakefile, sfhfile))
    if zinc is not None:
        command = '{0:s} -zinc'.format(command)
    if models is not None:
        command = '{0:s} -{1:s}'.format(command, models)
    if dAv is not None:
        command = '{0:s} -dAv={1:.2f}'.format(command, dAv)
    if dAvy is not None:
        command = '{0:s} -dAvy={1:.2f}'.format(command, dAvy)
    command = '{0:s} > {1:s}'.format(command, terminalfile)

    # Run calcsfh
    subprocess.call(command, shell=True)
    time.sleep(5)

    # Get fit value from calcsfh terminal output file
    with open(terminalfile, 'r') as f:
        for line in f:
            if line.startswith('Best fit: '):
                fit = float(line.split()[4].split('=')[1])

    # Get contents from .sfh and .sfh.cmd files
    with open(sfhfile, 'r') as f:
        sfh = f.read()
    with open(cmdfile, 'r') as f:
        sfhcmd = f.read()

    return fit, sfh, sfhcmd


def pattern_search(func, xy0, dxy, dxy0, divisor=2, xlim=None, ylim=None):
    """
    Find x and y that minimizes func using a 2-d pattern search.

    The pattern is a 9-point (3x3) array, starting with x0,y0 in the center
    and the eight surrounding points on a 2*dx by 2*dy square. func is
    evaluated at each point in the pattern and the point with the smallest
    result becomes the new central point. If the central point has the
    smallest result, then the step size is reduced; if the step size is
    already at the smallest allowed value, then the central point is the
    minimum of func.

    Parameters
    ----------
    func : function
        The function to be minimized.
    xy0 : tuple
        Initial seed point for the search as a pair of coordinates, (x, y).
    dxy : float or tuple
        Initial step size for point pattern. A single value is applied to
        both x and y; a tuple can be used for different x and y step sizes,
        (dx, dy).
    dxy0 : float or tuple
        Minimum value(s) for dxy.
    divisor : int or float, optional
        Amount by which to divide the step size (default is 2, i.e., divide
        step size in half).
    xlim, ylim : list or tuple, optional
        Two-element list containing the minimum and maximum values of x and
        y to search within (default is None). To set only one extreme of a
        coordinate and leave the other unconstrained, set the unconstrained
        extreme to +/-np.inf.

    Returns
    -------
    x, y, z :
        x and y coordinates where func is minimized, and the value of func
        at that point.

    """
    try:
        x0, y0 = xy0
    except TypeError:
        x0, y0 = xy0, xy0

    try:
        dx, dy = dxy
    except TypeError:
        dx, dy = dxy, dxy

    try:
        dx0, dy0 = dxy0
    except TypeError:
        dx0, dy0 = dxy0, dxy0

    if xlim is None:
        xlim = (-np.inf, np.inf)
    if ylim is None:
        ylim = (-np.inf, np.inf)

    while 1:
        # 9-point pattern
        x_list = np.array([x0, x0, x0+dx, x0+dx, x0+dx,
                           x0, x0-dx, x0-dx, x0-dx])
        y_list = np.array([y0, y0+dy, y0+dy, y0, y0-dy,
                           y0-dy, y0-dy, y0, y0+dy])

        # Keep x and y within limits
        i = ((xlim[0] <= x_list) & (x_list <= xlim[1]) &
             (ylim[0] <= y_list) & (y_list <= ylim[1]))
        x_list, y_list = x_list[i], y_list[i]

        # Evaluate func at each point in the pattern
        z_list = np.array([func(x, y) for x, y in zip(x_list, y_list)])

        # Find minimum point and either move or shrink the point pattern.
        # Break if central point is the minimum in the smallest pattern.
        i = z_list.argmin()
        if (x_list[i], y_list[i]) == (x0, y0):
            dx, dy = dx/divisor, dy/divisor
            if dx < dx0 and dy > dy0:
                dx = dx0
            elif dy < dy0 and dx > dx0:
                dy = dy0
            elif dx < dx0 and dy < dy0:
                break
        else:
            x0, y0 = x_list[i], y_list[i]

    return x_list[i], y_list[i], z_list[i]


def grid_search(func, xy1, xy2, nxy):
    """
    Find x and y that minimizes func using a 2-d grid search.

    The grid is specified by a pair of endpoints and a number of points, per
    coordinate (similar to numpy.linspace).

    Parameters
    ----------
    func : function
        The function to be minimized.
    xy1 : tuple
        Minimum x and y values of the grid as a pair of coordinates, (x, y).
    xy2 : tuple
        Maximum x and y values of the grid as a pair of coordinates, (x, y).
    nxy : tuple
        The number of x and y values in the grid (i.e., columns and rows)
        as a pair of integers, (nx, ny).

    Returns
    -------
    x, y, z :
        x and y coordinates where func is minimized, and the value of func
        at that point.

    """
    x1, y1 = xy1
    x2, y2 = xy2
    nx, ny = nxy

    # Set up the grid
    x, y = np.linspace(x1, x2, nx), np.linspace(y1, y2, ny)
    x_list, y_list = np.meshgrid(x, y)
    x_list, y_list = x_list.ravel(), y_list.ravel()

    # Evaluate func at each point in the grid
    z_list = np.array([func(x, y) for x, y in zip(x_list, y_list)])

    # Find minimum point in the grid
    i = z_list.argmin()

    return x_list[i], y_list[i], z_list[i]


def pointlist_search(func, xy_list):
    """
    Find x and y that minimizes func from a list of points.

    Parameters
    ----------
    func : function
        The function to be minimized.
    xy_list : list-like
        A list of points at which to evaluate func, one element per point.
        Each element is a two-element list of x, y coordinates for the
        point.

    Returns
    -------
    x, y, z :
        x and y coordinates where func is minimized, and the value of func
        at that point.

    """
    # Evaluate func at each point in the list
    z_list = np.array([func(x, y) for x, y in xy_list])

    # Find minimum point
    i = z_list.argmin()

    return xy_list[i][0], xy_list[i][1], z_list[i]


###TEST_CODE_BEGIN
### Experimental search algorithms
def _metropolis_hastings_search_test1(func, xy0, sigma, nsteps, dxy0=None,
                               xlim=None, ylim=None):
    """
    Find x and y that minimizes func using the Metropolis-Hastings algorithm.

    The location of each point is selected from a normal distribution with
    variance sigma**2 centered on the location of the previous point,
    creating a random walk. sigma therefore characterizes the step size or
    jump distance between successive points. sigma must be optimized for
    func such that it is large enough for "rapid mixing" (the walker
    doesn't spend too much time in a local minimum) and yet small enough
    for a reasonable acceptance rate (the walker isn't spending too much
    time jumping to low-probability areas). The result is the minimum value
    of func found after nstep interations.

    Because this function assumes fit values from calcsfh in MATCH where
    fit ~ -ln(p), the acceptance ratio r = p(x)/p(x_previous) may also be
    expressed as ln(r) = fit_previous - fit. If fit < fit_previous, then
    the current point is better than the previous point (p(x) > p(x_previous))
    and ln(r) >= 0 (r >= 1).

    Parameters
    ----------
    func : function
        The function to be minimized.
    xy0 : tuple
        Initial seed point for the search as a pair of coordinates, (x, y).
    dxy0 : float or tuple
        Minimum value(s) for dxy. Also determines the nearest multiple of x
        and y to which each step is rounded.
    sigma : float
        Width of the distribution from which to draw new points.
    nsteps : int
        Number of iterations after which the search is terminated and the
        minimum value of func found is returned.
    xlim, ylim : list or tuple, optional
        Two-element list containing the minimum and maximum values of x and
        y to search within (default is None). To set only one extreme of a
        coordinate and leave the other unconstrained, set the unconstrained
        extreme to +/-np.inf.

    Returns
    -------
    x, y, z :
        x and y coordinates where func is minimized, and the value of func
        at that point.

    """
    try:
        x0, y0 = xy0
    except TypeError:
        x0, y0 = xy0, xy0

    try:
        dx0, dy0 = dxy0
    except TypeError:
        dx0, dy0 = dxy0, dxy0

    if xlim is None:
        xlim = (-np.inf, np.inf)
    if ylim is None:
        ylim = (-np.inf, np.inf)

    x_list, y_list, z_list = [], [], []
    z0 = np.inf
    x, y = x0, y0
    for n in range(nsteps):
        z = func(x, y)
        x_list.append(x)
        y_list.append(y)
        z_list.append(z)

        lnr = (z0 - z)/2
        if 0 <= lnr or np.random.rand() <= np.exp(lnr):
            x0, y0, z0 = x, y, z

        while 1:
            x = np.round(np.random.normal(x0, sigma)/dx0) * dx0
            y = np.round(np.random.normal(y0, sigma)/dy0) * dy0
            if ((x != x0 or y != y0) and
                xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]):
                break

    x_list, y_list, z_list = np.array(x_list), np.array(y_list), np.array(z_list)
    i = z_list.argmin()

    return x_list[i], y_list[i], z_list[i]


def _metropolis_hastings_search_test2(func, xy0, sigma, niter, dxy0=None,
                               xlim=None, ylim=None):
    """
    Find x and y that minimizes func using the Metropolis-Hastings algorithm.

    The location of each point is selected from a normal distribution with
    variance sigma**2 centered on the location of the previous point,
    creating a random walk. sigma therefore characterizes the step size or
    jump distance between successive points. sigma must be optimized for
    func such that it is large enough for "rapid mixing" (the walker
    doesn't spend too much time in a local minimum) and yet small enough
    for a reasonable acceptance rate (the walker isn't spending too much
    time jumping to low-probability areas). The result is the minimum value
    that has been the minimum for the past niter iterations.

    Because this function assumes fit values from calcsfh in MATCH where
    fit ~ -ln(p), the acceptance ratio r = p(x)/p(x_previous) may also be
    expressed as ln(r) = fit_previous - fit. If fit < fit_previous, then
    the current point is better than the previous point (p(x) > p(x_previous))
    and ln(r) >= 0 (r >= 1).

    Parameters
    ----------
    func : function
        The function to be minimized.
    xy0 : tuple
        Initial seed point for the search as a pair of coordinates, (x, y).
    dxy0 : float or tuple
        Minimum value(s) for dxy. Also determines the nearest multiple of x
        and y to which each step is rounded.
    sigma : float
        Width of the distribution from which to draw new points.
    xlim, ylim : list or tuple, optional
        Two-element list containing the minimum and maximum values of x and
        y to search within (default is None). To set only one extreme of a
        coordinate and leave the other unconstrained, set the unconstrained
        extreme to +/-np.inf.

    Returns
    -------
    x, y, z :
        x and y coordinates where func is minimized, and the value of func
        at that point.

    """
    try:
        x0, y0 = xy0
    except TypeError:
        x0, y0 = xy0, xy0

    try:
        dx0, dy0 = dxy0
    except TypeError:
        dx0, dy0 = dxy0, dxy0

    if xlim is None:
        xlim = (-np.inf, np.inf)
    if ylim is None:
        ylim = (-np.inf, np.inf)

    x_list, y_list, z_list = [], [], []
    z0 = np.inf
    x, y, zmin = x0, y0, z0
    while 1:
        z = func(x, y)
        x_list.append(x)
        y_list.append(y)
        z_list.append(z)

        lnr = (z0 - z)/2
        if 0 <= lnr or np.random.rand() <= np.exp(lnr):
            x0, y0, z0 = x, y, z

        if z < zmin:
            zmincount = 0
            zmin = z
        zmincount += 1
        if zmincount > niter:
            break

        while 1:
            x = np.round(np.random.normal(x0, sigma)/dx0) * dx0
            y = np.round(np.random.normal(y0, sigma)/dy0) * dy0
            if ((x != x0 or y != y0) and
                xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]):
                break

    x_list, y_list, z_list = np.array(x_list), np.array(y_list), np.array(z_list)
    i = z_list.argmin()

    return x_list[i], y_list[i], z_list[i]
###TEST_CODE_END


def main():
    """
    The behavior of calcsfh and how the results are logged and organized
    are controlled here. The main tasks are,

    - Specify paths to directories and files and create any directories
      that do not yet exist.
    - Define a function, which gets passed to a search algorithm or run
      manager, that wraps around calcsfh to perform logging, retrieve
      results from previous calculations, and organize calcsfh output.
    - Specify the calcsfh parameters and command line options.
    - Call the search algorithm or run manager.

    """
    def func(Av, dAv):
        """
        Wrapper for calcsfh to record and recall fit values (preventing
        redundant calculations) and write .sfh and .cmd files to
        archive_dir.

        """
        # Retrieve prior result or start a new run
        Avstr, dAvstr = '{0:.2f}'.format(Av), '{0:.2f}'.format(dAv)
        if (Avstr, dAvstr) in record:
            fit = float(record[(Avstr, dAvstr)])
        else:
            ###TEST_CODE_BEGIN
            ### Do not run calcsfh in test mode; get fit value from a
            ### lookup table instead.
            if TEST:
                try:
                    fit = float(record2[(Avstr, dAvstr)])
                except KeyError:
                    fit = 1e10
            else:
                calcsfhargs = [photfile, fakefile, IMF, dmod, (Av, Av), res, logZ,
                               BF, badfrac, CMD, tbins]

                calcsfhkwargs = {'paramfile': paramfile, 'sfhfile': sfhfile,
                                 'cmdfile': cmdfile, 'terminalfile': terminalfile,
                                 'zinc': zinc, 'dAv': dAv, 'dAvy': dAvy,
                                 'models': models, 'path': path}

                fit, sfh, sfhcmd = calcsfh(*calcsfhargs, **calcsfhkwargs)
            ###TEST_CODE_END

            # Log the result, write sfh and cmd files for archive
            fitstr = '{0:.6f}'.format(fit)
            record[(Avstr, dAvstr)] = fitstr
            r_xml = etree.SubElement(data_xml, 'R')
            for val in [Avstr, dAvstr, fitstr]:
                c_xml = etree.SubElement(r_xml, 'C')
                c_xml.text = val
            format_table(table_xml)
            tree.write(logfile, encoding='UTF-8', xml_declaration=True)

            ###TEST_CODE_BEGIN
            ### Do not write calcsfh output files in test mode.
            if TEST:
                pass
            else:
                sfhout = os.path.join(archive_dir,
                                      '{0:s}_{1:s}.sfh'.format(Avstr, dAvstr))
                with open(sfhout, 'w') as f:
                    f.write(sfh)

                cmdout = '{0:s}.cmd'.format(sfhout)
                with open(cmdout, 'w') as f:
                    f.write(sfhcmd)
            ###TEST_CODE_END

        return fit


    name = sys.argv[1]


    # Path constants
    proj_dir = '/Users/Jake/Research/PHAT/UV_Regions/project/'
    phot_dir = os.path.join(proj_dir, 'phot')
    fake_dir = os.path.join(proj_dir, 'fake')
    sfh_dir = os.path.join(proj_dir, 'sfh')
    name_dir = os.path.join(sfh_dir, name)
    archive_dir = os.path.join(name_dir, 'archive')


    # Files
    photfile = os.path.join(phot_dir, '{0:s}.gst.match'.format(name))
    fakefile = os.path.join(fake_dir,
                            '12056_M31-B15-Fxx-WFC_F475W_F814W_gst.matchfake')
    paramfile = os.path.join(archive_dir, '{0:s}_temp.par'.format(name))
    sfhfile = os.path.join(archive_dir, '{0:s}_temp.sfh'.format(name))
    cmdfile = '{0:s}.cmd'.format(sfhfile)
    terminalfile = os.path.join(archive_dir, '{0:s}_temp.log'.format(name))
    logfile = os.path.join(name_dir, '{0:s}.xml'.format(name))


    # Check that paths exists where files are to be written
    for filename in [paramfile, sfhfile, cmdfile, terminalfile, logfile]:
        safe_mkdir(os.path.dirname(filename))

    # Check if an xml log file exists and create one if not
    try:
        tree = etree.parse(logfile)
    except IOError:
        table_xml = etree.Element('TABLE')
        coldef_xml = etree.SubElement(table_xml, 'COLDEF')
        data_xml = etree.SubElement(table_xml, 'DATA')

        coldef_list = [('A_V', 'float'), ('dA_V', 'float'), ('fit', 'float')]
        for cname, dtype in coldef_list:
            c_xml = etree.SubElement(coldef_xml, 'C')
            name_xml = etree.SubElement(c_xml, 'NAME')
            dtype_xml = etree.SubElement(c_xml, 'DTYPE')
            name_xml.text, dtype_xml.text = cname, dtype

        format_table(table_xml)
        tree = etree.ElementTree(table_xml)
        tree.write(logfile, encoding='UTF-8', xml_declaration=True)

    # Get results from previous runs
    table_xml = tree.getroot()
    data_xml = table_xml[1]
    record = {(row[0].text, row[1].text): row[2].text for row in data_xml}

    ###TEST_CODE_BEGIN
    ### Load lookup table of fit values; used instead of running calcsfh in
    ### test mode.
    if TEST:
        surfacefile  = ('/Users/Jake/Dropbox/Research/PHAT/UV_Regions/project/sfh'
                        '{0:s}/{0:s}_surface.xml'.format(name))
        tree2 = etree.parse(surfacefile)
        data2_xml = tree2.getroot()[1]
        record2 = {(row[0].text, row[1].text): row[2].text for row in data2_xml}
    ###TEST_CODE_END


    # Setup calcsfh parameters
    path = '/usr/local/MATCH/match2.4/bin/calcsfh'
    IMF = 'Kroupa'
    dmod = (24.47, 24.47)
    res = 0.05
    logZ = (-2.3, 0.1, 0.1)
    BF = 0.35
    badfrac = (1e-6, 1e-6)
    zinc = (-2.3, -0.9, -1.4, 0.1)
    dAvy = 0
    models = 'PADUA_AGB'

    filter1 = {'name': 'WFC475W', 'min': 16.0, 'max': 27.0}
    filter2 = {'name': 'WFC814W', 'min': 15.0, 'max': 26.2}
    CMD = {'V': filter1, 'I': filter2,
           'Vres': 0.1, 'V-I': (-0.5, 5.0, 0.05), 'fake_sm': 5,
           'exclude_gates': [[(1.25, 27.0), (5.0, 27.0),
                              (5.0, 21.0), (1.25, 21.0)]]}

    t1 = np.linspace(6.60, 9.00, (9.00-6.60)/0.05+1)
    t2 = np.linspace(9.10, 10.10, (10.10-9.10)/0.10+1)
    t = np.r_[t1, t2]
    tbins = zip(t[:-1], t[1:])


    # Run calcsfh in a given opration mode
    a = sys.argv
    mode = a[2]
    if mode == 'pattern':
        Av0, dAv0 = float(a[3]), float(a[4])  # Initial seed point
        d_Av, d_dAv = float(a[5]), float(a[6])  # initial step size
        d_Av0, d_dAv0 = float(a[7]), float(a[8])  # Minimum step size
        Avlim, dAvlim = (0.0, 1.6), (0.0, 5.0)  # Constrain Av and dAv
        pattern_search(func, (Av0, dAv0), (d_Av, d_dAv), (d_Av0, d_dAv0),
                       xlim=Avlim, ylim=dAvlim)

    elif mode == 'grid':
        Av1, dAv1 = float(a[3]), float(a[4])  # Lower Av, dAv limits
        Av2, dAv2 = float(a[5]), float(a[6])  # Upper Av, dAv limits
        nAv, ndAv = int(a[7]), int(a[8])  # Number dAv "rows" and Av "columns"
        grid_search(func, (Av1, dAv1), (Av2, dAv2), (nAv, ndAv))

    elif mode == 'list':
        pointfile = a[3]
        with open(pointfile, 'r') as f:
            pointlist = [line.split() for line in f if not line.isspace()]
        pointlist = np.array(pointlist, 'float')
        pointlist_search(func, pointlist)

    ###TEST_CODE_BEGIN
    ### Experimental search algorithms
    elif mode == 'Metropolis-Hastings':
        ###Av0, dAv0 = float(a[3]), float(a[4])  # Initial seed point
        ###sigma, nsteps = float(a[5]), int(a[6])  # Step size, num. iterations
        ###d_Av0, d_dAv0 = float(a[7]), float(a[8])  # Minimum step size
        ###Avlim, dAvlim = (0.0, 1.5), (0.0, 2.5)  # Constrain Av and dAv
        ###metropolis_hastings_search(func, (Av0, dAv0), sigma, nsteps,
        ###                           dxy0=(d_Av0, d_dAv0),
        ###                           xlim=Avlim, ylim=dAvlim)
        Av0, dAv0 = float(a[3]), float(a[4])  # Initial seed point
        sigma, niter = float(a[5]), int(a[6])  # Step size, num. iterations
        d_Av0, d_dAv0 = float(a[7]), float(a[8])  # Minimum step size
        Avlim, dAvlim = (0.0, 1.5), (0.0, 2.5)  # Constrain Av and dAv
        metropolis_hastings_search(func, (Av0, dAv0), sigma, niter,
                                   dxy0=(d_Av0, d_dAv0),
                                   xlim=Avlim, ylim=dAvlim)
    ###TEST_CODE_END


if __name__ == '__main__':
    main()
