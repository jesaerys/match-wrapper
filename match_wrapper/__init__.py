"""

===============
`match_wrapper`
===============

Python wrapper for MATCH, the SFH code by Andrew Dolphin [1]_.

The documentation for this package does not cover descriptions of the
utilities included in MATCH nor the details of their usage. Refer instead
to the MATCH 2.5 README and the utility help strings for that information.

`match_wrapper` requires MATCH to be installed (available from Andy; tested
with version 2.5). Python's `subprocess` module is used to run each
utility, so if calcsfh, zcombine, etc. work from the command line, then
`match_wrapper` should work as expected. The following packages are
required:

- `astropy <http://www.astropy.org>`_
- `numpy <http://www.numpy.org>`_

.. [1] Dolphin, A. E. 2002, MNRAS, 332, 91


Classes
-------

======================= ==================================================
|CalcsfhParam|          Class for storing calcsfh parameter file
                        information.
|CalcsfhParamFormatter| Formatter to assist in writing calcsfh parameter
                        files.
|CMDParam|              Class for storing CMD information.
|FlagFormatter|         Formatter to assist in writing flags to command
                        strings.
======================= ==================================================


Functions
---------

====================== ====================================================
|calcsfh|              calcsfh wrapper.
|hybridMC|             hybridMC wrapper.
|open_cmdfile|         Load data from a calcsfh ".cmd" file.
|open_zcbfile|         Load a zcombine/zcmerge output file into an
                       `astropy.table.Table`.
|write_zcbfile|        Write a zcombine/zcmerge output file from an
                       `astropy.table.Table`.
|write_zcombine_param| Write a zcombine parameter file for a given set of
                       age bin edges.
|zcmerge|              zcmerge wrapper.
|zcombine|             zcombine wrapper.
====================== ====================================================


============
Module Index
============

- `match_wrapper.match`
- `match_wrapper.io`
- `match_wrapper.util`


.. references

.. |FlagFormatter| replace:: `~match_wrapper.match.FlagFormatter`
.. |calcsfh| replace:: `~match_wrapper.match.calcsfh`
.. |hybridMC| replace:: `~match_wrapper.match.hybridMC`
.. |zcmerge| replace:: `~match_wrapper.match.zcmerge`
.. |zcombine| replace:: `~match_wrapper.match.zcombine`

.. |CMDParam| replace:: `~match_wrapper.io.CMDParam`
.. |CalcsfhParam| replace:: `~match_wrapper.io.CalcsfhParam`
.. |CalcsfhParamFormatter| replace:: `~match_wrapper.io.CalcsfhParamFormatter`
.. |open_cmdfile| replace:: `~match_wrapper.io.open_cmdfile`
.. |open_zcbfile| replace:: `~match_wrapper.io.open_zcbfile`
.. |write_zcbfile| replace:: `~match_wrapper.io.write_zcbfile`
.. |write_zcombine_param| replace:: `~match_wrapper.io.write_zcombine_param`

"""
from .match import (FlagFormatter, calcsfh, hybridMC, zcmerge, zcombine)
from .io import (CMDParam, CalcsfhParam, CalcsfhParamFormatter)
#from .util import *
