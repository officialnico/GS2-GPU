---
title: Response matrices
---

[TOC]

Introduction
============

GS2 calculates response matrices to be used as a part of the implicit time advance
algorithm. For large problems the calculation of these matrices can become expensive in
both time and memory. Here we discuss various topics of relevance to these matries.

Condition number
================

When using [field_option = 'local'](namelists/index.html#fields_knobs-field_option) we
calculate the condition number of the response matrix. There are several ways to view this
condition number:

-   Set `GK_VERBOSITY >= 1` to write the condition number to screen.
-   When running with [dump_response = T](namelists/index.html#fields_knobs-dump_response)
    we write the condition number into the netcdf file storing the response matrix. View
    with `ncdump -v condition <response_file> | tail -n 2`

When the condition number is large we consider the matrix to be ill-conditioned. In this
case we also display the condition number in the `<runname>.error` file and will write to
screen a warning if the module level parameter `debug` is true.

Saving and restoring the response matrices
==========================================

As the response matrices can be relatively expensive to calculate it is desirable to be
able to reuse these across restarted runs etc. Both `field_option = 'local'` and
`field_option = 'implicit'` offer the option to dump and read response matrices too/from
file. To write the response matrices set
[dump_response](namelists/index.html#fields_knobs-dump_response) true and to read existing
response matrices set [read_response](namelists/index.html#fields_knobs-read_response).

If `read_response` is true but the expected files are not found then the code will
fallback to calculating the response matrices directly. As such, it is recommended that
users set both `dump_response = T` and `read_response = T` for non-linear simulations. One
file is written per response matrix and we may have up to around `nx*ny` response matrices
per time step. The [response_dir](namelists/index.html#fields_knobs-response_dir) input
can be used to indicate a directory in which to write/read response matrices to/from,
avoiding cluttering the run directory.

The response matrix files generated can be reused on any processor count and are
compatible between `field_option = 'implicit'` and `field_option = 'local'` choices.
