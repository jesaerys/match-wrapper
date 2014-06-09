"""

====================
`match_wrapper.util`
====================

"""
# function to process large hybridMC output into separate .sfh or .zcb files.

# function to process large hybridMC output file into SFR PDF vs. time? A
# file with this data would be much more portable than the full data set
# output from MATCH - one would only need to sample the PDF evolution to
# build up a set of HMC SFHs.

# search stuff?
#
# It might be useful to rewrite UV_Regions/project/sfh/runcaclsfh.py,
# putting as much of the logging and search functionality in this module
# as possible. The updated runcalcsfh.py could then be an example script of
# how to use the match-wrapper package.
#
# How much of the searching can be handled by scipy.optimize?
