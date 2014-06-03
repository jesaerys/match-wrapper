"""

===============
`match_wrapper`
===============

Python wrapper for MATCH, the SFH code by Andrew Dolphin [1]_.

`match_wrapper` depends on MATCH (tested with version 2.5). Python's
`subprocess` module is used to run each utility, so if calcsfh, zcombine,
etc. work from the command line, then `match_wrapper` should also work. The
following packages are required:

- `astropy <http://www.astropy.org>`_
- `numpy <http://www.numpy.org>`_

.. [1] Dolphin, A. E. 2002, MNRAS, 332, 91


Functions
---------

========== =================
|calcsfh|  calcsfh wrapper.
|zcmerge|  zcmerge wrapper.
|zcombine| zcombine wrapper.
========== =================


============
Module Index
============

- `match_wrapper.match`
- `match_wrapper.io`
- `match_wrapper.util`


.. references

.. |calcsfh| replace:: `~match_wrapper.match.calcsfh`
.. |zcmerge| replace:: `~match_wrapper.match.zcmerge`
.. |zcombine| replace:: `~match_wrapper.match.zcombine`

"""
from .match import (calcsfh, zcmerge, zcombine)
from .io import open
#from .util import *
