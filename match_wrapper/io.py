"""

==================
`match_wrapper.io`
==================

One liner.

Description.

classes
-------

======================= ==================================
`CMDParam`              Class
`CalcsfhParam`          Class
`CalcsfhParamFormatter` Class for storing CMD information.
======================= ==================================


functions
---------


"""
import numpy as np


class CMDParam(object):

    """Class for storing CMD information.

    Most parameter and attribute names are taken from the MATCH 2.5 README.
    The name from the MATCH 2.5 README is given in quotes if the parameter
    or attribute name is different.

    Each parameter initializes a corresponding attribute and defaults to
    None unless stated otherwise.

    Parameters
    ----------
    Vname : str, optional
    Iname : str, optional
    Vstep : float, optional
    VImin : float, optional
    VImax : float, optional
    VIstep : float, optional
    fake_sm : int, optional
    exclude_gates : list, optional
        Default is an empty list.
    combine_gates : list, optional
        Default is an empty list.

    Attributes
    ----------
    Vname : str
        "V". `CalcsfhParam.filters` must have a filter dictionary with the
        same name.
    Iname : str
        "I". `CalcsfhParam.filters` must have a filter dictionary with the
        same name.
    Vstep : float
    VImin : float
        "V-Imin"
    VImax : float
        "V-Imax"
    VIstep : float
        "V-Istep"
    fake_sm : float
    exclude_gates : list
        List of exclude gates. Each gate is a set of coordinates of the
        form ``[(x1, y1), (x2, y2), (x3, y3), (x4, y4)]``.
    combine_gates : list
        List of combine gates (see `exclude_gates`).

    Methods
    -------
    Nexclude_gates
        Property (get only). Length of the `exclude_gates` list.
    Ncombine_gates
        Property (get only). Length of the `combine_gates` list.

    """

    def __init__(self, **kwargs):
        self.Vname = kwargs.get('Vname')
        self.Iname = kwargs.get('Iname')
        self.Vstep = kwargs.get('Vstep')
        self.VImin = kwargs.get('VImin')
        self.VImax = kwargs.get('VImax')
        self.VIstep = kwargs.get('VIstep')
        self.fake_sm = kwargs.get('fake_sm')

        gates = kwargs.get('exclude_gates')
        try:
            gates[0][0][0]
        except TypeError:
            try:
                gates[0][0]
            except TypeError:
                gates = []
            else:
                gates = [gates]
        self.exclude_gates = gates

        gates = kwargs.get('combine_gates')
        try:
            gates[0][0][0]
        except TypeError:
            try:
                gates[0][0]
            except TypeError:
                gates = []
            else:
                gates = [gates]
        self.combine_gates = gates

    @property
    def Nexclude_gates(self):
        return len(self.exclude_gates)

    @property
    def Ncombine_gates(self):
        return len(self.combine_gates)


class CalcsfhParam(object):

    """Class for storing calcsfh parameter file information.

    Most parameter and attribute names are taken from the MATCH 2.5 README.
    The name from the MATCH 2.5 README is given in quotes if the parameter
    or attribute name is different.

    Each parameter initializes a corresponding attribute and defaults to
    None unless stated otherwise.

    Parameters
    ----------
    IMF : {None, 'Kroupa', 'Salpeter', float}, optional
        'Kroupa' is equivalent to -1.0, 'Salpeter' is equivalent to 1.35.
    dmodmin : float, optional
    dmodmax : float, optional
    Avmin : float, optional
    Avmax : float, optional
    step : float, optional
    logZmin : float, optional
        See `mode`.
    logZmax : float, optional
        See `mode`.
    logZstep : float, optional
        See `mode`.
    logZimin : float, optional
        See `mode`.
    logZimax : float, optional
        See `mode`.
    logZfmin : float, optional
        See `mode`.
    logZfmax : float, optional
        See `mode`.
    logZspread : float, optional
        See `mode`.
    BF :  float, optional
    Bad0 : float, optional
    Bad1 : float, optional
    CMDs : CMDParam or list, optional
        A `CMDParam` instance, or a list of one or more `CMDParam`
        instances. Default is an empty list.
    filters : list, optional
        A list of dictionaries for the filters referenced by the `CMDParam`
        instances in `CMDs`. Default is an empty list.
    agebins : list, optional
        Default is an empty list.
    linage : bool, optional
        Default is False.
    logZcentral : list, optional
        Default is a list of length `Ntbins` where each element is None.
    SFR : list, optional
        Default is a list of length `Ntbins` where each element is None.
    bgCMDs : list, optional
        A list of dictionaries, one for each background/foreground CMD.
        Default is an empty list.
    mode : {None, 'zinc', 'setz'}, optional

    Attributes
    ----------
    IMF : float
    dmodmin : float
        "m-Mmin".
    dmodmax : float
        "m-Mmax".
    Avmin : float
    Avmax : float
    step : float
        Both "d(m-M)" and "dAv".
    logZmin : float
        See `mode`.
    logZmax : float
        See `mode`.
    logZstep : float
        "dlogZ"; see `mode`.
    logZimin : float
        Minimum initial (oldest age) metallicity for 'zinc' mode; see `mode`.
    logZimax : float
        Maximum initial (oldest age) metallicity for 'zinc'; see `mode`.
    logZfmin : float
        Minimum final (youngest age) metallicity for 'zinc'; see `mode`.
    logZfmax : float
        Maximum final (youngest age) metallicity for 'zinc'; see `mode`.
    logZspread : float
        Metallicity spread for 'setz' mode; see `mode`.
    BF : float
    Bad0 : float
    Bad1 : float
    CMDs : list
        List of `CMDParam` instances, one per CMD.
    filters : list
        List of filters referenced by the `CMDParam` instances in `CMDs`.
        Each filter is a dictionary containing the keys,

        - 'name': "V" or "I"
        - 'min': Bright magnitude limit; "Vmin" or "Imin"
        - 'max': Faint magnitude limit; "Vmax" or "Imax"

        A filter corresponds to V or I depending on whether it is the bluer
        or redder filter in a given CMD.
    agebins : list
        List of edges of the age bins (either in yr or as log10(t/yr),
        depending on `linage`). The ith and i+1th elements correspond to
        "To" (youngest edge) and "Tf" (oldest edge) of the ith age bin.
    linage : bool
        True if the `agebins` values are in linear years, False if the
        values are ``log10(age/yr)``.
    logZcentral : list
        Central metallicity values of the age bins for 'setz' mode; see
        `mode`. The ith value corresponds to the ith age bin, and a value
        must be specified for each age bin. Length is `Ntbins`.
    SFR : list
        Force the SFR in particular age bins to the values in this list.
        The ith element corresponds to the ith age bin, and, if not None,
        appears in the parameter file as the last number in the line for
        the age bin. SFR is only fixed where the value in the list is not
        None. Length is `Ntbins`.

        .. note:: This feature is not documented in the MATCH 2.5 README!

    bgCMDs : list
        List of background/foreground CMD dictionaries. Each dictionary
        contains the keys,

        - 'nbins': Size of the smoothing kernel.
        - 'scale': Set the scaling for the background CMD. If negative, a
          variable number of stars is used.
        - 'filename': Optional; path to the file containing the background
          CMD data. If None (default), a smoothed version of the observed
          CMD is used.
        - 'cmdfile': Optional; True if 'filename' is formatted like a
          ".cmd" file output by calcsfh. Default is False, i.e., 'filename'
          has two-columns like an input photometry file for calcsfh.

    mode : str
        The following attributes are either required or ignored depending
        on the mode:

        None (default):

        - required: `logZmin`, `logZmax`, `logZstep`
        - ignored: `logZimin`, `logZimax`, `logZfmin`, `logZfmax`,
          `logZspread`, `logZcentral`

        'zinc':

        - required: `logZmin`, `logZmax`, `logZstep`, `logZimin`,
          `logZimax`, `logZfmin`, `logZfmax`
        - ignored: `logZspread`, `logZcentral`

        'setz':

        - required: `logZspread`, `logZcentral`
        - ignored: `logZmin`, `logZmax`, `logZstep`, `logZimin`,
          `logZimax`, `logZfmin`, `logZfmax`

    Methods
    -------
    Ncmds
        Property (get only). Length of the `CMDs` list.
    To
        Property (get only).
    Tf
        Property (get only).
    Ntbins
        Property (get only). One less than the length of the `agebins` list
        (length of `To` and `Tf`).
    read
        Create a `CalcsfhParam` instance from a calcsfh parameter file.
    write
        Write to a calcsfh parameter file.

    """

    def __init__(self, **kwargs):
        IMF = kwargs.get('IMF')
        if IMF == 'Kroupa':
            IMF = -1.0
        elif IMF == 'Salpeter':
            IMF = 1.35
        self.IMF = IMF

        self.dmodmin = kwargs.get('dmodmin')
        self.dmodmax = kwargs.get('dmodmax')
        self.Avmin = kwargs.get('Avmin')
        self.Avmax = kwargs.get('Avmax')
        self.step = kwargs.get('step')

        self.logZmin = kwargs.get('logZmin')  # default or zinc; not setz
        self.logZmax = kwargs.get('logZmax')  # default or zinc; not setz
        self.logZstep = kwargs.get('logZstep')  # default or zinc; not setz
        self.logZimin = kwargs.get('logZimin')  # zinc only
        self.logZimax = kwargs.get('logZimax')  # zinc only
        self.logZfmin = kwargs.get('logZfmin')  # zinc only
        self.logZfmax = kwargs.get('logZfmax')  # zinc only
        self.logZspread = kwargs.get('logZspread')  # setz only

        self.BF = kwargs.get('BF')
        self.Bad0 = kwargs.get('Bad0')
        self.Bad1 = kwargs.get('Bad1')

        CMDs = kwargs.get('CMDs', [])
        self.CMDs = [CMDs] if isinstance(CMDs, CMDParam) else CMDs
        self.filters = kwargs.get('filters', [])

        self.agebins = kwargs.get('agebins', [])
        self.linage = kwargs.get('linage', False)
        self.logZcentral = kwargs.get('logZcentral', [None]*self.Ntbins)
        self.SFR = kwargs.get('SFR', [None]*self.Ntbins)

        bgCMDs = kwargs.get('bgCMDs', [])
        try:
            bgCMDs['nbins']  # just a test for a dict or a list of dict
        except TypeError:
            pass
        else:
            bgCMDs = [bgCMDs]
        self.bgCMDs = bgCMDs

        self.mode = kwargs.get('mode')

    @property
    def Ncmds(self):
        return len(self.CMDs)

    @property
    def To(self):
        return self.agebins[:-1]

    @property
    def Tf(self):
        return self.agebins[1:]

    @property
    def Ntbins(self):
        l = len(self.agebins)-1
        return l if l>0 else 0

    def read(self, filename):
        """Create a `CalcsfhParam` instance from a calcsfh parameter file.

        Parameters
        ----------
        filename : str
            Absolute path to the input parameter file.

        Returns
        -------
        CalcsfhParam

        """
        with open(filename, 'r') as f:
            lines = f.readlines()

        items = []
        n = 0

        # IMF, dmod, Av parameters
        keys = ['IMF', 'dmodmin', 'dmodmax', 'step', 'Avmin', 'Avmax']
        vals = [float(val) for val in lines[n].split()[:-1]]
        items += zip(keys, vals)
        n += 1

        # Metallicity parameters
        vals = [float(val) for val in lines[n].split()]
        if len(vals) == 1:
            mode = 'setz'
            keys = ['logZspread']
        elif len(vals) == 7:
            mode = 'zinc'
            keys = ['logZmin', 'logZmax', 'logZstep', 'logZimin',
                    'logZimax', 'logZfmin', 'logZfmax']
        else:
            mode = None
            keys = ['logZmin', 'logZmax', 'logZstep']
        items += zip(keys, vals)
        n += 1

        # Binary fraction and upper/lower bad fractions
        keys = ['BF', 'Bad0', 'Bad1']
        vals = [float(val) for val in lines[n].split()]
        items += zip(keys, vals)
        n += 1

        # CMDs (part 1 of 2)
        Ncmds = int(lines[n])
        n += 1
        CMDs, filternames = [], []
        for i in range(Ncmds):
            vals = lines[n].split()
            Vname, Iname = vals[5].split(',')
            filternames += [Vname, Iname]
            CMDitems = [
                ('Vstep', float(vals[0])),
                ('VIstep', float(vals[1])),
                ('fake_sm', int(vals[2])),
                ('VImin', float(vals[3])),
                ('VImax', float(vals[4])),
                ('Vname', Vname),
                ('Iname', Iname)
                ]
            CMDs.append(CMDitems)
            n += 1

        # Filter list
        filters = []
        while lines[n].split()[-1] in filternames:
            vals = lines[n].split()
            fdict = {'min': float(vals[0]), 'max': float(vals[1]), 'name': vals[2]}
            filters.append(fdict)
            n += 1
        items.append(('filters', filters))

        # Gate list
        for i in range(Ncmds):
            vals = lines[n].split()

            Nexc = int(vals[0])
            if Nexc > 0:
                points = vals[1:Nexc*8+1]
                gates = []
                for j in range(Nexc):
                    xy = [float(pnt) for pnt in points[j*8:(j+1)*8]]
                    gates.append(zip(xy[0::2], xy[1::2]))
                CMDs[i].append(('exclude_gates', gates))

            vals = vals[Nexc*8+1]
            Ncom = int(vals[0])
            if Ncom > 0:
                points = vals[1:Ncom*8+1]
                gates = []
                for j in range(Nexc):
                    xy = [float(pnt) for pnt in points[j*8:(j+1)*8]]
                    gates.append(zip(xy[0::2], xy[1::2]))
                CMDs[i].append(('combine_gates', gates))

            n += 1

        # CMDs (part 2 of 2)
        CMDs = [CMDParam(**dict(CMDitems)) for CMDitems in CMDs]
        items.append(('CMDs', CMDs))

        # Age bins
        Ntbins = int(lines[n])
        n += 1
        bins = [tuple(line.split()) for line in lines[n:n+Ntbins]]
        age1, age2, logZcentral, SFR = [], [], [], []
        sfrcol = 3 if mode == 'setz' else 2
        for vals in bins:
            vals = [float(val) for val in vals]
            age1.append(vals[0])
            age2.append(vals[1])
            if mode == 'setz':
                logZcentral.append(vals[2])
            else:
                logZcentral.append(None)
            if len(vals) == sfrcol+1:
                SFR.append(vals[sfrcol])
            else:
                SFR.append(None)
        agebins = age1 + [age2[-1]]
        linage = True if agebins[0] < 0 else False
        keys = ['agebins', 'linage', 'logZcentral', 'SFR']
        vals = [agebins, linage, logZcentral, SFR]
        items += zip(keys, vals)
        n += Ntbins

        # Background/foreground CMDs
        Nbgcmds = len(lines) - n
        if Nbgcmds > 0:
            bgCMDs = []
            for i in range(Nbgcmds):
                vals = lines[n].split()
                bgdict = {}
                nbins = int(vals[1])
                if nbins < 0:
                    bgdict['cmdfile'] = True
                    nbins *= -1
                bgdict['nbins'] = nbins
                bgdict['scale'] = float(vals[2]) if '.' in vals[2] else int(vals[2])
                if len(vals) == 4:
                    bgdict['filename'] = vals[3]
                bgCMDs.append(bgdict)
            items.append(('bgCMDs', bgCMDs))

        items.append(('mode', mode))

        return CalcsfhParam(**dict(items))

    def write(self, filename, formatter=None):
        """Write to a calcsfh parameter file.

        Parameters
        ----------
        filename : str
            Absolute path to the output parameter file.
        formatter : CalcsfhParamFormatter or function, optional
            Any function that takes a parameter name (`key`) and a value
            (`val`) as the first and second arguments, and returns a string
            representation of the value. `CalcsfhParamFormatter` is used by
            default.

        Returns
        -------
        None

        """
        if formatter is None:
            formatter = CalcsfhParamFormatter()

        # IMF, dmod, Av parameters
        pars = ['IMF', 'dmodmin', 'dmodmax', 'step', 'Avmin', 'Avmax', 'step']
        line = [formatter(key, self.__getattribute__(key)) for key in pars]
        line = ' '.join(line)
        lines = [line]

        # Metallicity parameters
        if self.mode == 'zinc':
            pars = ['logZmin', 'logZmax', 'logZstep',
                    'logZimin', 'logZimax', 'logZfmin', 'logZfmax']
        elif self.mode == 'setz':
            pars = ['logZspread']
        else:
            pars = ['logZmin', 'logZmax', 'logZstep']
        line = [formatter(key, self.__getattribute__(key)) for key in pars]
        line = ' '.join(line)
        lines.append(line)

        # Binary fraction and upper/lower bad fractions
        pars = ['BF', 'Bad0', 'Bad1']
        line = [formatter(key, self.__getattribute__(key)) for key in pars]
        line = ' '.join(line)
        lines.append(line)

        # Number of CMDs
        line = formatter('Ncmds', self.Ncmds)
        lines.append(line)

        # CMD list
        pars = ['Vstep', 'VIstep', 'fake_sm', 'VImin', 'VImax', 'Vname', 'Iname']
        for CMD in self.CMDs:
            line = [formatter(key, CMD.__getattribute__(key)) for key in pars]
            line = '{0:s} {1:s},{2:s}'.format(' '.join(line[:-2]), line[-2], line[-1])
            lines.append(line)

        # Filter list
        pars = ['min', 'max', 'name']
        for filt in self.filters:
            line = [formatter(key, filt[key]) for key in pars]
            line = ' '.join(line)
            lines.append(line)

        # Gate list
        for CMD in self.CMDs:
            Nexc = formatter('Nexclude_gates', CMD.Nexclude_gates)
            exc = ' '.join([formatter('exclude_gates', x)
                            for gate in CMD.exclude_gates
                            for point in gate for x in point])
            line1 = '{0:s} {1:s}'.format(Nexc, exc) if exc else Nexc
            Ncom = formatter('Ncombine_gates', CMD.Ncombine_gates)
            com = ' '.join([formatter('combine_gates', x)
                            for gate in CMD.combine_gates
                            for point in gate for x in point])
            line2 = '{0:s} {1:s}'.format(Ncom, com) if com else Ncom
            line = '{0:s} {1:s}'.format(line1, line2)
            lines.append(line)

        # Number of age bins
        line = formatter('Ntbins', self.Ntbins)
        lines.append(line)

        # Age bins
        linage = -1 if self.linage else 1
        for i in range(self.Ntbins):
            edge1 = formatter('To', linage*self.To[i])
            edge2 = formatter('Tf', linage*self.Tf[i])
            logZc = formatter('logZcentral', self.logZcentral[i]) if self.mode == 'setz' else ''
            SFR = formatter('SFR', self.SFR[i]) if self.SFR[i] else ''
            row = [edge1, edge2, logZc, SFR]
            line = ' '.join(val for val in row if val)  # non-empty strs only
            lines.append(line)

        # Background/foreground CMDs
        for CMD in self.bgCMDs:
            cmdfile = -1 if CMD.get('cmdfile') else 1
            nbins = formatter('nbins', cmdfile*CMD['nbins'])
            scale = formatter('scale', CMD['scale'])
            fname = CMD.get('filename', '')
            row = ['-1', nbins, scale, fname]
            line = ' '.join(val for val in row if val)  # non-empty strs only
            # Space, or no space, between scale and filename?
            lines.append(line)

        with __builtins__.open(filename, 'w') as f:
            f.writelines('{:s}\n'.format(line) for line in lines)

        return None


class CalcsfhParamFormatter(object):

    """Formatter to assist in writing calcsfh parameter files.

    The `CalcsfhParam.write` method formats each value in the output
    parameter file using a formatter function that takes a key and a value
    and returns a string. `CalcsfhParamFormatter` has a call method that
    looks up a format string based on the key, and uses it to format the
    value. Each key has a corresponding attribute in `CalcsfhParam`. The
    format strings may be adjusted from their default values, shown in the
    Parameters section below.

    Parameters
    ----------
    IMF : '{ :.2f}', optional
    dmodmin : '{ :.2f}', optional
    dmodmax : '{ :.2f}', optional
    Avmin : '{ :.2f}', optional
    Avmax : '{ :.2f}', optional
    step : '{ :.2f}', optioanl
    logZmin : '{ :.1f}', optional
    logZmax : '{ :.1f}', optional
    logZstep : '{ :.1f}', optional
    logZimin : '{ :.1f}', optional
    logZimax : '{ :.1f}', optional
    logZfmin : '{ :.1f}', optional
    logZfmax : '{ :.1f}', optional
    logZspread : '{ :.1f}', optional
    BF : '{ :.2f}', optional
    Bad0 : '{ :.6f}', optional
    Bad1 : '{ :.6f}', optional
    Ncmds : '{ :d}', optional
    Vstep : '{ :.2f}', optional
    VIstep : '{ :.2f}', optional
    fake_sm : '{ :d}', optional
    VImin : '{ :.2f}', optional
    VImax : '{ :.2f}', optional
    Vname : '{ :s}', optional
    Iname : '{ :s}', optional
    min : '{ :.2f}', optional
    max : '{ :.2f}', optional
    name : '{ :s}', optional
    Nexclude_gates : '{ :d}', optional
    exclude_gates : '{ :.2f}', optional
    Ncombine_gates : '{ :d}', optional
    combine_gates : '{ :.2f}', optional
    Ntbins : '{ :d}', optional
    To : '{ :.2f}', optional
    Tf : '{ :.2f}', optional
    logZcentral : '{ :.1f}', optional
    SFR : '{ :.3e}', optional
    nbins : '{ :d}', optional
    scale : '{ :d}', optional
    filename : '{ :s}', optional

    Attributes
    ----------
    fmt_dict : dict
        Dictionary of format strings for all values in a calcsfh parameter
        file.

    Methods
    -------
    __call__
        Return a string for the given key and value.

    """

    def __init__(self, **kwargs):
        fmt_dict = {
            'IMF': '{:.2f}',
            'dmodmin': '{:.2f}',
            'dmodmax': '{:.2f}',
            'Avmin': '{:.2f}',
            'Avmax': '{:.2f}',
            'step': '{:.2f}',
            'logZmin': '{:.1f}',
            'logZmax': '{:.1f}',
            'logZstep': '{:.1f}',
            'logZimin': '{:.1f}',
            'logZimax': '{:.1f}',
            'logZfmin': '{:.1f}',
            'logZfmax': '{:.1f}',
            'logZspread': '{:.1f}',
            'BF': '{:.2f}',
            'Bad0': '{:.6f}',
            'Bad1': '{:.6f}',
            'Ncmds': '{:d}',
            'Vstep': '{:.2f}',
            'VIstep': '{:.2f}',
            'fake_sm': '{:d}',
            'VImin': '{:.2f}',
            'VImax': '{:.2f}',
            'Vname': '{:s}',
            'Iname': '{:s}',
            'min': '{:.2f}',
            'max': '{:.2f}',
            'name': '{:s}',
            'Nexclude_gates': '{:d}',
            'exclude_gates': '{:.2f}',
            'Ncombine_gates': '{:d}',
            'combine_gates': '{:.2f}',
            'Ntbins': '{:d}',
            'To': '{:.2f}',
            'Tf': '{:.2f}',
            'logZcentral': '{:.1f}',
            'SFR': '{:.3e}',
            'nbins': '{:d}',
            'scale': '{:d}',
            'filename': '{:s}',
            }
        for key, val in kwargs.items():
            if key in fmt_dict:
                fmt_dict[key] = val

    def __call__(self, key, val):
        """Return a string for the given key and value.

        Parameters
        ----------
        key : str
            The key corresponding to a format string in `fmt_dict`.
        val :
            The value to be formatted.

        Returns
        -------
        str
            The formatted value.

        """
        return self.fmt_dict[key].format(val)


def _parse_sfhfile(filename):
    return None


def _parse_cmdfile(filename):
    return None


def _parse_zcbfile(filename):
    return None


def open(filename, kind=None):
    return None





# read/write zcombine parameter file


def test():
    V = dict(name='WFC475W', min=16.00, max=27.00)
    I = dict(name='WFC814W', min=15.00, max=26.20)
    filters = [V, I]

    gate_x = [1.25, 5.00, 5.00, 1.25]
    gate_y = [27.00, 27.00, 21.00, 21.00]
    CMD = CMDParam(Vname='WFC475W', Iname='WFC814W', Vstep=0.10, VImin=-0.50,
                   VImax=5.00, VIstep=0.05, fake_sm=5,
                   exclude_gates=zip(gate_x, gate_y))

    logages1 = np.linspace(6.60, 9.00, (9.00-6.60)/0.05+1)
    logages2 = np.linspace(9.10, 10.10, (10.10-9.10)/0.10+1)
    agebins = np.hstack((logages1, logages2))

    param = CalcsfhParam(
            IMF='Salpeter',
            dmodmin=24.47, dmodmax=24.47, Avmin=0.45, Avmax=0.45, step=0.05,
            logZmin=-2.3, logZmax=0.1, logZstep=0.1,
            logZimin=-2.3, logZimax=-0.9, logZfmin=-1.4, logZfmax=0.1,
            BF=0.35, Bad0=1e-6, Bad1=1e-6,
            CMDs=CMD, filters=filters, agebins=agebins, mode='zinc')

    param.write('/Users/Jake/Desktop/testpar.par')


if __name__ == '__main__':
    pass
    #test()


