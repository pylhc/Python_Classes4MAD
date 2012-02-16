#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      05 Sep 2006
#$Revision: 1.1 $
"""
    Plotting imports for PyDSTool, from Matplotlib.

    Robert Clewley, March 2006.
"""

# --------------------------------------------------------------------------
# VERSION HISTORY
#
# v1.1, August 2006
#   Added save_fig function to save figure in multiple formats.
#
# --------------------------------------------------------------------------

__all__ = ['show', 'plot', 'pylab', 'save_fig']

import matplotlib
ver = matplotlib.__version__.split(".")
try:
    if int(ver[0]) == 0 and int(ver[1]) < 65:
        import matplotlib.matlab as pylab
    else:
        import matplotlib.pylab as pylab
except RuntimeError, err:
    if str(err) == 'could not open display':
        failed=True
    else:
        raise
else:
    failed=False

if failed:
    # Dummy plot overrides for PyDSTool when matplotlib fails to import
    def show():
        print "Warning: show does not work!"

    def plot(*args, **kw):
        print "Warning: plot does not work!"

    def save_fig(fignum, fname, formats=[]):
        print "Warning: plot does not work!"

    print "Warning: matplotlib failed to import properly and so is not"
    print "  providing a graphing interface"
    pylab = None   # will cause an error if someone tries to access in order to plot
else:
    import os
    from Trajectory import Trajectory

    # Mini-utility to solve annoying pylab vs. IDLE problem for displaying windows
    def show():
        if os.name == 'nt' or os.name == 'dos' or os.name == 'ce':
            pylab.draw()
        else:
            pylab.show()

    # Convenient shorthand to permit singleton numeric types and Trajectories
    # in the plot arguments without first converting them to lists or arrays.
    def plot(*args, **kw):
        new_args = list(args)
        if type(args[0]) in [int, float]:
            new_args[0] = [args[0]]
        elif type(args[0]) == Trajectory:
            try:
                new_args[0] = args[0].sample()
            except:
                raise RuntimeError, "Could not sample trajectory with default options for plotting"
        if len(args) > 1:
            if type(args[1]) in [int, float]:
                new_args[1] = [args[1]]
            elif type(args[1]) == Trajectory:
                try:
                    new_args[1] = args[1].sample()
                except:
                    raise RuntimeError, "Could not sample trajectory with default options for plotting"
        return pylab.plot(*tuple(new_args), **kw)
    
    def save_fig(fignum, fname, formats=['png','svg','eps']):
        """Save figure fignum to multiple files with different formats
        and extensions given by the formats argument.
        These are platform-dependent and are specific to matplotlib's support.
        """
        for f in formats:
            pylab.figure(fignum).savefig(fname+'.'+f)
