"""
.. module: plotting

This module contains easy to use macros for presentable plots using gnuplot. A simple example usage is::

    from plotting import BetaBeat
    dir1='path_to_my_old_results'
    dir2='path_to_my_new_results'
    BetaBeat(dir1,dir2)

.. moduleauthor:: Yngve Inntjore Levinsen <yngve.inntjore.levinsen@cern.ch>
"""

from BetaBeat import BetaBeat
from Dispersion import Dispersion
from NormDispersion import NormDispersion
from Coupling import Coupling
from _append import set_append
