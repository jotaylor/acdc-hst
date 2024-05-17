..
   Note: Items in this toctree form the top-level navigation. See `api.rst` for the `autosummary` directive, and for why `api.rst` isn't called directly.

.. toctree::
   :hidden:

   Home <self>
   What is ACDC? <about>
   Installation and Setup <install>
   Databases <databases>
   Usage <usage>
   Package Index <_autosummary/acdc>

Another COS Background Correction (ACDC)
========================================

**Why another COS background correction?**

COS spectroscopic science in the extreme UV regime is limited by detector background noise. The COS pipeline, calcos does not perform an optimized background subtraction. In order to achieve the maximum scientific value of the COS instrument, we have a designed a custom characterization and correction of the COS FUV dark rate.

With this package we can:

* create and maintain databases needed to measure the dark rate as a function of time, position, PHA, and more
* create COS/FUV superdarks
* use superdarks to perform custom dark corrections
* analyze the efficacy of custom dark-corrected COS data

