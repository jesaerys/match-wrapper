"""

==================
`match_wrapper.io`
==================

"""
import numpy as np

from .util import Param


class CMDParam(object):

    """Class for storing CMD information.

    The name of a parameter in the MATCH 2.5 README is quoted if it differs
    from the corresponding property or attribute name below.

    Parameters
    ----------
    Vname : str, optional
        Initialize the `Vname` attribute. Default is None.
    Iname : str, optional
        Initialize the `Iname` attribute. Default is None.
    Vstep : float, optional
        Initialize the `Vstep` attribute. Default is None.
    VImin : float, optional
        Initialize the `VImin` attribute. Default is None.
    VImax : float, optional
        Initialize the `VImax` attribute. Default is None.
    VIstep : float, optional
        Initialize the `VIstep` attribute. Default is None.
    fake_sm : int, optional
        Initialize the `fake_sm` attribute. Default is None.
    exclude_gates : list, optional
        Initialize the `exclude_gates` attribute. Default is an empty list.
    combine_gates : list, optional
        Initialize the `combine_gates` attribute. Default is an empty list.

    Attributes
    ----------
    Vname : str
        "V"; `CalcsfhParam.filters` must have a filter dictionary with the
        same name.
    Iname : str
        "I"; `CalcsfhParam.filters` must have a filter dictionary with the
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

    The name of a parameter in the MATCH 2.5 README is quoted if it differs
    from the corresponding property or attribute name below.

    Parameters
    ----------
    IMF : {None, 'Kroupa', 'Salpeter', float}, optional
        Initialize the `IMF` attribute. 'Kroupa' is equivalent to -1.0,
        'Salpeter' is equivalent to 1.35. Default is None.
    dmodmin : float, optional
        Initialize the `dmodmin` attribute. Default is None.
    dmodmax : float, optional
        Initialize the `dmodmax` attribute. Default is None.
    Avmin : float, optional
        Initialize the `Avmin` attribute. Default is None.
    Avmax : float, optional
        Initialize the `Avmax` attribute. Default is None.
    step : float, optional
        Initialize the `step` attribute. Default is None.
    logZmin : float, optional
        Initialize the `logZmin` attribute; see `mode`. Default is None.
    logZmax : float, optional
        Initialize the `logZmax` attribute; see `mode`. Default is None.
    logZstep : float, optional
        Initialize the `logZstep` attribute; see `mode`. Default is None.
    logZimin : float, optional
        Initialize the `logZimin` attribute; see `mode`. Default is None.
    logZimax : float, optional
        Initialize the `logZimax` attribute; see `mode`. Default is None.
    logZfmin : float, optional
        Initialize the `logZfmin` attribute; see `mode`. Default is None.
    logZfmax : float, optional
        Initialize the `logZfmax` attribute; see `mode`. Default is None.
    logZspread : float, optional
        Initialize the `logZspread` attribtue; see `mode`. Default is None.
    BF :  float, optional
        Initialize the `BF` attribute. Default is None.
    Bad0 : float, optional
        Initialize the `Bad0` attribute. Default is None.
    Bad1 : float, optional
        Initialize the `Bad1` attribute. Default is None.
    CMDs : CMDParam or list, optional
        Initialize the `CMDs` attribute. The list should contain one or
        more `CMDParam` instances. A single `CMDParam` instance is
        converted to a length-1 list. Default is None.
    filters : list, optional
        Initialize the `filters` attribute. The list should contain
        dictionaries for the filters referenced by the CMDs in `CMDs`.
        Default is an empty list.
    agebins : list, optional
        Initialize the `agebins` attribute. Default is an empty list.
    linage : bool, optional
        Initialize the `linage` attribute. Default is False.
    logZcentral : list, optional
        Initialize the `logZcentral` attribute. Default is a list of length
        ``len(agebins)-1`` where each element is None.
    SFR : list, optional
        Initialize the `SFR` attribute. Default is a list of length
        ``len(agebins)-1`` where each element is None.
    bgCMDs : list, optional
        Initialize the `bgCMDs` attribute. The list should contain one
        dictionary for each background CMD. Default is an empty list.
    mode : {None, 'zinc', 'setz'}, optional
        Initialize the `mode` attribute. Default is None.

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
        Minimum initial metallicity for 'zinc' mode; see `mode`.
    logZimax : float
        Maximum initial metallicity for 'zinc'; see `mode`.
    logZfmin : float
        Minimum final metallicity for 'zinc'; see `mode`.
    logZfmax : float
        Maximum final metallicity for 'zinc'; see `mode`.
    logZspread : float
        Metallicity spread for 'setz' mode; see `mode`.
    BF : float
    Bad0 : float
    Bad1 : float
    CMDs : list
        List of `CMDParam` instances, one per CMD.
    filters : list
        List of filters referenced by the CMDs in `CMDs`. Each filter is a
        dictionary containing the keys,

        - 'name': "V" or "I"
        - 'min': Bright magnitude limit; "Vmin" or "Imin"
        - 'max': Faint magnitude limit; "Vmax" or "Imax"

        A filter corresponds to V or I depending on whether it is the bluer
        or redder filter in a given CMD.
    agebins : list
        List of edges of the age bins (either in yr or as log10(t/yr),
        depending on `linage`). The ith and i+1th elements correspond to
        "To" and "Tf" of the ith age bin.
    linage : bool
        If True, values in `agebins` are linear years. Default is False
        (log10).
    logZcentral : list
        Central metallicity values of the age bins for 'setz' mode; see
        `mode`. The ith value corresponds to the ith age bin, and a value
        must be specified for each age bin. Length is ``len(agebins) - 1``.
    SFR : list
        Force the SFR in particular age bins to the values in this list.
        The ith element corresponds to the ith age bin, and, if not None,
        appears in the parameter file as the last number in the line for
        the age bin. SFR is only fixed where the value in the list is not
        None.

        .. note:: This feature is not documented in the MATCH 2.5 README!

    bgcmd : list
        List of background CMD dictionaries. Each dictionary contains the
        keys,

        - 'nbins': Size of the smoothing kernel
        - 'scale': Set the scaling for the background CMD. If negative, a
          variable number of stars is used.
        - 'filename': Optional; path to the file containing the background
          CMD data. If None (default), a smoothed version of the observed
          CMD is used.
        - 'cmdfile': Optional; True if 'filename' is formatted like a .cmd
          file output by calcsfh. Default is False, i.e., 'filename' is
          like a two-column input photometry file for calcsfh.

    mode : str
        The following attributes are required or ignored depending on the mode:

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
        Return a `CalcsfhParam` instance from a calcsfh parameter file.
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
            bgCMDs['nbins']
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
        """
        # Defaults
        linage = False
        agebins_col3 = False
        bgcmd = None

        with open(filename, 'r') as f:
            lines = f.readlines()

        line = lines[0].split()

        return CalcsfhParam(IMF, dmod, Av, step, logZ, BF, badfrac, CMD,
                            agebins, linage=linage,
                            agebins_col3=agebins_col3, bgcmd=bgcmd)
        """
        return None

    def write(self, filename, formatter=None):
        """Write to a calcsfh parameter file.

        Parameters
        ----------
        filename : str
            Absolute path to the output parameter file.
        formatter : function, optional
            Any function that takes a parameter name (`key`) and a value
            (`val`) as the first and second arguments, and returns a string
            representation of the value. If None (default), then the
            returned value is ``str(val)``.

        Returns
        -------
        None

        """
        def default_formatter(key, val):
            return str(val)

        if formatter is None:
            formatter = default_formatter

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
            edge1 = formatter('To', linage*self.agebins[i])
            edge2 = formatter('Tf', linage*self.agebins[i+1])
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
            filename = CMD.get('filename', '')
            row = ['-1', nbins, scale, filename]
            line = ' '.join(val for val in row if val)  # non-empty strs only
            # Space, or no space, between scale and filename?
            lines.append(line)

        with open(filename, 'w') as f:
            f.writelines('{:s}\n'.format(line) for line in lines)

        return None


def checkparam(obj):
    return obj if isinstance(obj, Param) else Param(obj)


def _parse_sfhfile(filename):
    return None


def _parse_cmdfile(filename):
    return None


def _parse_zcbfile(filename):
    return None


def open(filename, kind=None):
    return None





# read/write zcombine parameter file

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


def custom_formatter(key, val):
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
            'Bad0': None,
            'Bad1': None,
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
    if key in ['Bad0', 'Bad1']:
        n = abs(int(np.log10(val)))  # the expoonent is negative
        result = '{1:{0}f}'.format(n, val)
    else:
        result = fmt_dict[key].format(val)

    return result




fmt = '{:.2f}'
(cp.CMD[0].V.min, cp.CMD[0].V.max,
cp.CMD[0].I.min, cp.CMD[0].I.max,
cp.CMD[0].Vstep,
cp.CMD[0].cmin, cp.CMD[0].cmax, cp.CMD[0].cstep,
cp.CMD[0].exclude_gates[0],
cp.IMF, cp.dmod, cp.Av, cp.step, cp.logZ, cp.BF, cp.badfrac,
cp.tbins)





