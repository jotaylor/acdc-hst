Usage
=====

The ACDC package has many uses- custom dark correction, data analysis,
database creation and querying, and superdark creation. Some common use
examples are provided below.

Custom Dark Correction
----------------------

To run the custom dark correction, you must have the ACDC superdarks.
TODO- document where to obtain these files.

Now, simply run by executing the following command:

::

  acdc -i <indir> -o <outdir>

Where ``<indir>`` is the input directory with corrtags you wish to calibrate,
and ``<outdir>`` is the directory where custom-calibrated products should be written.

Data Analysis
-------------

TODO

Superdark Creation
------------------

You must have the SQLite database downloaded to create superdarks.
TODO- document where DB can be downloaded from.

To create the pre-defined quiescent and active superdarks for each
segment+HV combination, simply run:

::

  make_all_superdarks

If you wish to create a custom superdark, you will need to do so interactively
in a ipython/python session or create a script with the following commands:

::

  from acdc.superdark.cos_fuv_superdark import Superdark
  S = Superdark(hv=163, segment="FUVA", mjdstarts=[58000], mjdends=[58220], 
                bin_pha=1, phastart=1, phaend=31, outfile="my_custom_superdark")

You may fill in the arguments as you see fit.

Database Querying
-----------------

Many common queries have been pre-written for convenience. 
Some examples for queries on the SQLite database are below.

**Return all dark events**

Returns expstart, solar_flux, latitude, longitude, segment, hv, region, saa_distance,
and all dark_pha<x> columns

::

  from acdc.database.query_cos_dark import all_darks
  df = all_darks()

**Return filenames based on MJD, HV, and segment**

Returns distinct filenames

::

  from acdc.database.query_cos_dark import files_by_mjd
  df = files_by_mjd(mjdstart, mjdend, segment, hv)

**Return number of dark events based on MJD**

Returns expstart, solar_flux, latitude, longitude, segment, hv, region, saa_distance,
and all dark_pha<x> columns

::

  from acdc.database.query_cos_dark import counts_by_mjd
  df = counts_by_mjd(mjdstart, mjdend)

Some examples for queries on the MySQL database are below.

**Return all dark events for a single HV**

Unless otherwise specified, returns all columns. 
I.e. ``returncols`` is 
``[xcorr, ycorr, pha, mjd, hv, segment, filename, and proposid]``.

::

  from acdc.database.query_dark_events import all_rows
  df = all_rows(hv, returncols)

**Execute a SQL query**

Execute a SQL command string.

::

  from acdc.database.query_dark_events import sql_query
  sql_query(hv, sqlquery)
  # E.g. sql_query(163, "SELECT * FROM darkeventshv163 WHERE id=1")

**Select rows based on column values**

This is the equivalent of a ``WHERE col=val`` SQL command.
Unless otherwise specified, returns all columns. 
I.e. ``returncols`` is 
``[xcorr, ycorr, pha, mjd, hv, segment, filename, and proposid]``.

::

  from acdc.database.query_dark_events import equals_query
  equals_query(hv, returncols, **kwargs)


**Select rows based on a range**

Possible range to query across are: MJD, X position, Y position, and PHA.

::

  from acdc.database.query_dark_events import range_query
  range_query(163, mjdstart=58484, mjdend=58490)
  range_query(163, x0=1000, x1=5000, y0=200, y1=512)
  range_query(163, pha0=20, pha1=25)
