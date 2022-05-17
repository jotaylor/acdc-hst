# Another COS Background Correction (ACDC) ⚡

## Another background correction?
COS spectroscopic science in the extreme UV regime is limited by detector 
background noise. The COS pipeline, [`calcos`](https://github.com/spacetelescope/calcos) 
does not perform an optimized 
background subtraction. In order to achieve the maximum scientific value of
the COS instrument, we have a designed a custom characterization and correction
of the COS FUV dark rate.

With this package we can: 
* create and maintain databases needed to measure the dark rate as a function of time, position, PHA, and more
* create COS/FUV superdarks
* use superdarks to perform custom dark corrections
* analyze the efficacy of custom dark-corrected COS data

## Installation
If you do not already have Conda installed, you need to download and install
either Miniconda or Anaconda. Miniconda provides a bare minimum Conda
environment. Anaconda provides a full Conda root environment along with
many other tools, libraries, and utilities.
* get [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
* get [Anaconda](https://www.anaconda.com/products/individual)

Now clone this repository on your local machine. 

For SSH connection: `git clone git@github.com:spacetelescope/cos_dark.git`

For HTTPS connection: `git clone https://github.com/spacetelescope/cos_dark.git`

#### Installing into a fresh environment
`cd` into the cloned repository and execute:

```
conda create -n <env_name> python=<version>
conda activate <env_name>
pip install .
```

where `<env_name>` is the name of the environment that will be created.
You need at least python version 3.9, so fill in `<version>` with whatever
version above >=3.9 that you want.
This only installs dependencies required by the `acdc` packae. You will
need to install anything else manually, e.g.:

```
pip install ipython multiprocessing
```

#### Installing into an existing environment
Activate your desired environment, then `cd` into the cloned repository and install using `pip`:

```
pip install ".[docs]"
```

The extra `[docs]` ensures that you pick up the requirements needed to build the docs locally.

## Usage

For full usage instructions, refer to the documentation. While the repository is private,
the docs must be built locally. To do so, navigate into the cloned repository and into
the `docs/` directory. Once there, run the following command:

```
make html
```

It will take a minute or two to build all the docs. Once finished, you can open the 
docs in your default web browser with the following command:

```
open _build/html/index.html 
```

From there you can click and navigate the webpage as if it were hosted online normally.
