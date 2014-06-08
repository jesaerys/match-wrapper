"""

=====================
`match_wrapper.match`
=====================

Wrappers for MATCH utilities.

Refer to the MATCH 2.5 README for full documentation on each command. Like
their command line counterparts, the functions in this module read and
write data in files, rather than processing data arrays stored in memory.


Classes
-------

=============== ========================================================
`FlagFormatter` Formatter to assist in writing flags to command strings.
=============== ========================================================


Functions
---------

========== =================
`calcsfh`  calcsfh wrapper.
`hybridMC` hybridMC wrapper.
`zcmerge`  zcmerge wrapper.
`zcombine` zcombine wrapper.
========== =================

"""
import subprocess
import time


class FlagFormatter(object):

    """Formatter to assist in writing flags to command strings.

    Each MATCH utility wrapper formats a flag for the command string using
    a formatter function that takes a key and value and returns a string.
    `FlagFormatter` has a call method that looks up a format string based
    on the key, and uses it to format the value. Each key corresponds to an
    argument in a wrapper function. The format strings may be adjusted from
    their default values, shown in the Parameters section below.

    A key can be assigned a list or tuple of format strings to specify
    individual formats for each value in a multivalued flag (e.g.,
    `diskAv`); otherwise the same format is used for all values.

    Parameters
    ----------

    Attributes
    ----------
    fmt_dict : dict
        Dictionary of format strings for flags used by MATCH utilities.

    Methods
    -------
    __call__
        Return a string for the given key and value.

    """

    def __init__(self, **kwargs):
        self.fmt_dict = {
            'dAv': '{:.2f}',
            'diskAv': '{:.2f}',
            'dAvy': '{:.2f}',
            'zerobin': '{:s}',
            'incmult': '{:.2f}',
            'incmultsig': '{:.2f}',
            'errmult': '{:.2f}',
            'errmultsig': '{:.2f}',
            'logterr': '{:.2f}',
            'logterrsig': '{:.2f}',
            'mbolerr': '{:.2f}',
            'mbolerrsig': '{:.2f}',
            'alpha': '{:.2f}',
            }
        for key, val in kwargs.items():
            if key in self.fmt_dict:
                self.fmt_dict[key] = val

    def __call__(self, key, val, delim=','):
        """Return a string for the given key and value.

        Parameters
        ----------
        key : str
            The name of the value, and also the key corresponding to a
            format string in `fmt_dict`.
        val :
            The value to be formatted. No formatting is actually applied if
            `val` is True, False, or None.
        delim : str, optional
            Delimiter for multivalued flags. Default is ','.

        Returns
        -------
        str or None
            None is returned if `val` is False or None. If `val` is True,
            then it is ignored and the returned string is of the form
            '-key'. Otherwise, `val` is formatted based on `key` and the
            returned string is of the form '-key=val'. If `val` is
            multivalued, the subvalues are separatd by `delim`.

        """
        if not val:  # ignore None and False
            result = None
        elif val is True:  # True is a special case
            result = '-{:s}'.format(key)
        else:
            fmt = self.fmt_dict[key]
            if self._islistlike(val):
                if self._islistlike(fmt):
                    valstr_list = [f.format(v) for v, f in zip(val, fmt)]
                else:
                    valstr_list = [fmt.format(v) for v in val]
                valstr = delim.join(valstr_list)
            else:
                valstr = fmt.format(val)
            result = '-{0:s}={1:s}'.format(key, valstr)

        return result

    def _islistlike(self, obj):
        """True if the object is iterable like a list and is *not* a string."""
        return ((hasattr(obj, '__iter__') or hasattr(obj, '__getitem__')) and
                not isinstance(obj, basestring))


def calcsfh(paramfile, photfile, fakefile, sfhfile, **kwargs):
    """calcsfh wrapper.

    Most command line flags in the MATCH 2.5 README have a corresponding
    keyword argument in the call signature. Mutually exclusive options,
    such as certain mode flags and flags for model and transformation
    selections, are available as values for the `mode`, `model`, and
    `transformation` keywords.

    The keyword arguments default to None, or False if boolean, unless
    stated otherwise.

    Parameters
    ----------
    paramfile : str
        Path to the calcsfh parameter file.
    photfile : str
        Path to the photometry file.
    fakefile : str
        Path to the fake photometry file.
    sfhfile : str
        Path to the SFH output file. The best-fit model CMD data is
        automatically written to the same path with a ".cmd" suffix.
    outfile : str, optional
        Path to a file to capture stdout.
    mode : {None, 'ddist', 'setz', 'ssp', 'zinc'}, optional
        Set the mode in which calcsfh is run (there are other modes that
        can also be specified, but only one of this particular set can be
        used at a time).
    full : bool, optional
    verb : bool, optional
    allstars : bool, optional
    mcdata : bool, optional
    dAv : float, 3-tuple, optional
    diskAv : 4-tuple, 7-tuple, optional
    dAvy : float, optional
    zerobin : str, optional

        .. note:: Partially supported. I'm not sure what the correct form
           of the argument is, so it is currently limited to a preformatted
           string. The value of the zerobin flag in the command string will
           show up exactly as it is given here.

    incmult : float, optional
    incmultsig : float, optional
    errmult : float, optional
    errmultsig : float, optional
    logterr : float, optional
    logterrsig : float, optional
    mbolerr : float, optional
    mbolerrsig : float, optional
    model : {None, 'BASTInov', 'BASTI', 'BASTI02nov', 'BASTI02', 'DART',
            'PADUA00', 'PADUA06', 'VICTORIA', 'VICTORIASS'}, optional
        Override the default stellar evolution model. Default is Padua 2006
        with C/O ratio support.
    alpha : {-0.2, 0, 0.2, 0.2}, optional
    transformation : {None, 'BASTInov', 'BASTI', 'BASTI02nov', 'BASTI02',
            'DART', 'PADUA00', 'PADUA06', 'VICTORIA', 'VICTORIASS'},
            optional
        Override the default transformation for the selected model. Default
        is to use the "new" Girardi transformations.
    formatter : FlagFormatter or function, optional
        Any function that takes a flag name (`key`) and a value (`val`) as
        the first and second arguments, and returns a string representation
        of the value. Formatting for multivalued flags (e.g., `diskAv`)
        need to be handled as well. `FlagFormatter` is used by default.
    norun : bool, optional
        Do not actually run calcsfh if True, just return the command
        string. This is useful for checking the command before running it.

    Returns
    -------
    str or None
        The calcsh command string is returned if `norun` is True, otherwise
        the return value is None.

    """
    cmd = ('calcsfh {0:s} {1:s} {2:s} {3:s}'
           .format(paramfile, photfile, fakefile, sfhfile))

    flag_list = []

    # Mode flag
    modeflag = kwargs.get('mode')
    if modeflag:
        kwargs[modeflag] = True
        flag_list.append(modeflag)

    # Other mode flags
    flag_list += ['full', 'verb', 'allstars', 'mcdata']

    # Model generation flags
    flag_list += ['dAv', 'diskAv', 'dAvy', 'zerobin', 'incmult',
                  'incmultsig', 'errmult', 'errmultsig', 'logterr',
                  'logterrsig', 'mbolerr', 'mbolerrsig']

    # Model selection
    modelflag = kwargs.get('model')
    if modelflag:
        kwargs[modelflag] = True
        flag_list.append(modelflag)
        if modelflag in ['DART']:  # `alpha` for Dartmouth models only
            flag_list.append('alpha')

    # Transformation selection
    transflag = kwargs.get('transformation')
    if transflag:
        kwargs[transflag] = True
        flag_list.append(transflag)

    # Flag strings
    formatter = kwargs.get('formatter', FlagFormatter())
    flagstr_gen = (formatter(flag, kwargs.get(flag)) for flag in flag_list)
    flagstr_list = [flagstr for flagstr in flagstr_gen if flagstr]  # filter out None
    if flagstr_list:
        cmd = '{0:s} {1:s}'.format(cmd, ' '.join(flagstr_list))

    # stdout file
    outfile = kwargs.get('outfile')
    if outfile:
        cmd = '{0:s} > {1:s}'.format(cmd, outfile)

    norun = kwargs.get('norun')
    if norun:
        result = cmd
    else:
        # Run calcsfh
        subprocess.call(cmd, shell=True)
        time.sleep(5)  # Finish writing potentially large output files
        result = None
    return result


def hybridMC(inputfile, outputfile, **kwargs):
    """hybridMC wrapper.

    ***Most command line flags in the MATCH 2.5 README have a corresponding
    keyword argument in the call signature. Mutually exclusive options,
    such as certain mode flags and flags for model and transformation
    selections, are available as values for the `mode`, `model`, and
    `transformation` keywords.

    ***The keyword arguments default to None, or False if boolean, unless
    stated otherwise.

    Parameters
    ----------
    formatter : FlagFormatter or function, optional
        Any function that takes a flag name (`key`) and a value (`val`) as
        the first and second arguments, and returns a string representation
        of the value. ***Formatting for multivalued flags (e.g., `diskAv`)
        need to be handled as well. `FlagFormatter` is used by default.
    norun : bool, optional
        Do not actually run hybridMC if True, just return the command
        string. This is useful for checking the command before running it.

    Returns
    -------
    str or None
        The hybridMC command string is returned if `norun` is True,
        otherwise the return value is None.

    """
    return None


def zcombine():
    return None


def zcmerge():
    return None

#combine1
#fake
#makefake
#sspcombine
#stats

if __name__ == '__main__':
    pass
