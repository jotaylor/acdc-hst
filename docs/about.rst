What is ACDC?
=======================

Do we really need **A**\ nother **C**\ OS **D**\ ark **C**\ orrection?

COS spectroscopic science in the extreme UV regime is limited by detector background noise. The COS pipeline, calcos does not perform an optimized background subtraction. In order to achieve the maximum scientific value of the COS instrument, we have a designed a custom characterization and correction of the COS FUV dark rate.

With this package we can:

* create and maintain databases needed to measure the dark rate as a function of time, position, PHA, and more
* create COS/FUV superdarks
* use superdarks to perform custom dark corrections
* analyze the efficacy of custom dark-corrected COS data

Custom Superdarks
-----------------

For each valid inflight segment and HV combination, at least two custom superdarks
exist. One is a superdark made from individual dark exposures obtained during a 
"quiescent" period (sometimes referred to as "low period"). 
The other is made from dark exposures obtained during an
"active" period (sometimes referred to as "high period"). 
The time limits for quiescent and active superdarks for each segment+HV combination
were manually selected based on analysis using data from the 
``cos_dark`` database (see :ref:`sqlite_db`).


Custom Dark Correction
----------------------

The custom dark correction takes as input a directory of COS/FUV corrtag files. 
Based on the observational configuration of these files, appropriate custom
superdarks are automatically selected. For example, if the input corrtag files
are all FUVA data of High Voltage 163, then the FUVA HV=163 quiescent and 
active superdarks will be used during calibration.

Next, for each input corrtag, an exposure-specific model superdark is created
on-the-fly. This model, or predicted, superdark is created by comparing
the non-PSA and non-WCA regions of each corrtag to the quiescent and active
superdarks. Each superdark is scaled and then combined to create a model
superdark for each corrtag. 

The custom dark correction is then performed directly on the corrtags.
For each event in the corrtag, the coordinates and PHA value are
checked to determine if it is within the limits used to construct the 
superdarks. By default, these spatial limits are defined as the active
area of the detector and the PHA limits are defined as >=2 and <=23. 
If an event is within these limits, then the modeled dark current at the
corresponding location is subtracted from the science event's ``EPSILON``
value. This serves as the custom dark correction.

Finally, custom dark-corrected corrtags are then run through CalCOS to
produce 1D spectra. However, to avoid performing a dark correction twice,
the ``BACKCORR`` switch must be set to ``OMIT``. Unfortunately, the 
twozone extraction algorithm does not work with ``BACKCORR=OMIT``. As a result,
``XTRACTALG`` must be set to use the boxcar extraction algorithm. Consequently,
all twozone switches must be set to ``OMIT`` as well (``ALGNCORR``, ``TRCECORR``).
The final result is a 1D spectrum that has been extracted using the boxcar
method. 


The Finer Details
-----------------

For an in-depth explanation of the COS dark behavior, creation of superdarks,
application of the custom dark correction, and more, please see
:download:`our full technical report <COS_Background_Full_Report.pdf>`.
