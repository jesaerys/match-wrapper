match-wrapper
=============

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

Documentation build instructions (requires `sphinx <http://sphinx-doc.org/>`_
and `numpydoc <https://github.com/numpy/numpydoc>`_)::

  cd match-wrapper/docs
  make html
  open _build/html/index.html

.. [1] Dolphin, A. E. 2002, MNRAS, 332, 91
