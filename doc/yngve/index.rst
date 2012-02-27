

Yngve's library
***************

This module contains some simple code which you may or may not find useful.



Function yngve.tfs
------------------

The tfs module can be used to load a tfs table into memory. It is created as 
two lookup dictionaries, one for the summary table and one for the main table.
The main table is the first return argument.

Usage example::

    from yngve import tfs
    t_table,sum_table=tfs('my_twiss.tfs')

Automatic documentation
=======================

.. autofunction:: yngve.tfs


Module yngve.mad2placet
-----------------------

This submodule holds functions for converting Mad-X output to Placet input.

Usage example::

    from yngve import tfs,mad2placet

    t_table=tfs('my_errors.tfs')[0]
    mad2placet.errors(t_table)

Automatic documentation
=======================

.. automodule:: yngve.mad2placet
    :members:
