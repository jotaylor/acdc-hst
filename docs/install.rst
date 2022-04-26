Installation and Setup
======================

If you do not already have Conda installed, you need to download and install
either Miniconda or Anaconda. Miniconda provides a bare minimum Conda
environment. Anaconda provides a full Conda root environment along with
many other tools, libraries, and utilities.

* get `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
* get `Anaconda <https://www.anaconda.com/products/individual>`_

Now clone this repository on your local machine. 

For SSH connection: ``git clone git@github.com:spacetelescope/cos_dark.git``

For HTTPS connection: ``git clone https://github.com/spacetelescope/cos_dark.git``

Installing into a fresh environment
-----------------------------------
``cd`` into the cloned repository, and execute:

::

  conda create -n <env_name> python
  conda activate <env_name>
  pip install .

where ``<env_name>`` is the name of the environment that will be created.
This only installs dependencies required by the ``acdc`` package. You will
need to install anything else manually, e.g.:

::

  pip install ipython multiprocessing

Installing into an existing environment
---------------------------------------

Activate your desired environment, then `cd` into the cloned repository and install using ``pip``:

::

  pip install .

.. _env_vars:

Environment Variables
---------------------

You must define two environment variables that are required by the code. The first is 
is the directory that houses COS reference files:

::

  export lref=<path>

where ``<path>`` is the path to your local CRDS cache.

The second is the directory that houses the custom COS/FUV superdarks. If you have not
already downloaded the superdarks, refer to :ref:`download_darks`. Once downloaded, the 
environment variable is defined as:

::

  export ACDC_SUPERDARKS=<path>

where ``<path>`` is the superdark directory.

If defined on the command line, these variables must be defined in every terminal.
You may want to consider adding the definitions to your local .bashrc, .bash_profile,
or similar file.
