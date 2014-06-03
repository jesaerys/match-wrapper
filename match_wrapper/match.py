"""

=====================
`match_wrapper.match`
=====================

Wrappers for MATCH utilities.

Refer to the MATCH 2.5 README for full documentation on each command. Like
their command line counterparts, the functions in this module read and
write data in files, rather than processing data arrays stored in memory.

Functions
---------

========== =================
`calcsfh'  calcsfh wrapper.
`zcmerge'  zcmerge wrapper.
`zcombine' zcombine wrapper.
========== =================

"""
import subprocess
import time

from . import Param


def _make_flagstr(flag, val):
    """Make a formatted string for the given flag and value."""
    param = Param(val)

    if not param.val:  # ignore None and False
        flagstr = None
    elif param.val is True:  # True is a special case
        flagstr = '-{:s}'.format(flag)
    else:
        flagstr = '-{0:s}={1:s}'.format(flag, param)

    return flagstr


def calcsfh(paramfile, photfile, fakefile, sfhfile, **kwargs):
    """calcsfh wrapper.

    Most command line flags have a corresponding keyword argument. Mutually
    exclusive options, such as certain mode flags and flags for model and
    transformation selection, are available as values for the `mode`,
    `model`, and `transformation` keywords. Any keyword value (including an
    individual value in a list) can have an optional format string to
    control how it is formatted in the calcsfh command. A format string is
    paired with its corresponding value in a tuple, e.g., ::

      # unspecified formats
      calcsfh(..., dAv=0.1)
      calcsfh(..., dAv=(0, 0.67, 1))

      # specified formats
      calcsfh(..., dAv=(0.1, '{:.2f}'))
      calcsfh(..., dAv=[0, (0.67, '{:.2f}'), 1])  # format 2nd value only

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
        used at a time). Default is None.
    full : bool, optional
        Default is False.
    verb : bool, optional
        Default is False.
    allstars : bool, optional
        Default is False.
    mcdata : bool, optional
        Default is False.
    dAv : float, 3-tuple, optional
        Default is None.
    diskAv : 4-tuple, 7-tuple, optional
        Default is None.
    dAvy : float, optional
        Default is None.
    zerobin : str, optional
        Default is None.

        .. note:: Partially supported. I'm not sure what the correct form
           of the argument is, so it is currently limited to a preformatted
           string. The value of the zerobin flag in the command string will
           show up exactly as it is given here.

    incmult : float, optional
        Default is None.
    incmultsig : float, optional
        Default is None.
    errmult : float, optional
        Default is None.
    errmultsig : float, optional
        Default is None.
    logterr : float, optional
        Default is None.
    logterrsig : float, optional
        Default is None.
    mbolerr : float, optional
        Default is None.
    mbolerrsig : float, optional
        Default is None.
    model : {None, 'BASTInov', 'BASTI', 'BASTI02nov', 'BASTI02', 'DART',
            'PADUA00', 'PADUA06', 'VICTORIA', 'VICTORIASS'}, optional
        Override the default stellar evolution model. Default is None
        (Padua 2006 with C/O ratio support).
    alpha : {-0.2, 0, 0.2, 0.2}, optional
        Default is None.
    transormation : {None, 'BASTInov', 'BASTI', 'BASTI02nov', 'BASTI02',
            'DART', 'PADUA00', 'PADUA06', 'VICTORIA', 'VICTORIASS'},
            optional
        Override the default transformation for the selected model. Default
        is None (the "new" Girardi transformations).
    norun : bool, optional
        Do not actually run calcsfh, just return the command. This is
        useful for checking the command before running it.

    Returns
    -------
    str or None
        If `norun` is True, the calcsh command string is returned.

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
    flagstr_gen = (_make_flagstr(flag, kwargs.get(flag)) for flag in flag_list)
    flagstr_list = [flagstr for flagstr in flagstr_gen if flagstr]
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


def zcombine():
    return None


def zcmerge():
    return None


if __name__ == '__main__':
    pass
