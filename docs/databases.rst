Databases
=========

Currently, two databases are maintained in order to characterize the 
FUV dark rate.

.. _sqlite_db:

The SQLite Database, ``cos_dark``
---------------------------------

This database
tracks the dark counts over binned regions of the FUV detectors. 
The SQLite database is available on Box- contact Jo Taylor for details.

To access the database, you will need 
`SQLite <https://www.sqlite.org/index.html>`_
if you do not already have it. Once you have it, simply type:

::

  sqlite3 cos_dark.db

and start exploring! For tips on SQLite queries, click 
`here <https://www.tutorialspoint.com/sqlite/sqlite_select_query.htm>`_.

Interactive plots of dark counts vs. time are available in the 
`plot_dark_vs_time.ipynb` jupyter notebook.

SQLite Database Format
~~~~~~~~~~~~~~~~~~~~~~

The format of the ``cos_darks`` SQLite database is as follows:

.. csv-table:: SQLite Darks Table
   :file: sqlite_darks.csv
   :widths: 25, 25, 25
   :header-rows: 1

.. csv-table:: SQLite Solar Table
   :file: sqlite_solar.csv
   :widths: 25, 25, 25
   :header-rows: 1

Solar flux values are obtained from `NOAA <https://www.swpc.noaa.gov/phenomena/f107-cm-radio-emissions>`_.

The MySQL database, ``hstcal``
------------------------------

The MySQL database, ``hstcal``, is hosted on an
STScI internal server. To access the ``hstcal`` database, you need MySQL version >= 8. 
For connection information, contact Jo Taylor.

MySQL Database Format
~~~~~~~~~~~~~~~~~~~~~

There is a table for each HV setting, and each table is named like
``darkeventshv<hvsettings>``. ``<hvsetting>`` is one of the following:

* 163
* 167
* 169
* 171
* 173
* 178
* other

The format of all ``hstcal`` tables is:

.. csv-table:: MySQL Darkevents Tables
   :file: mysql_darkeventshv.csv
   :widths: 25, 25, 25
   :header-rows: 1

