.. _system_stacking:

********************************************
Making System Stacking Fingerprints (SSFs)
********************************************

.. currentmodule:: stacker.pairwise_distance

Identifying a Pi-Stacking Pair
------------------------------
Users can identify a residue pair engaging in pi-stacking
by using the :ref:`System Stacking Fingerprint Pipeline <system_stacking>`.

It is recommended that you parallelize for large trajectories, it calculates the SSF
per frame in parallel.

Comparing Pi-Stacking Pairs
----------------------------