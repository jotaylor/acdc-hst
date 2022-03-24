# Another COS Background Correction (ACDC)

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
`cd` into the cloned repository, and execute:

```
conda create -n <env_name> python
conda activate <env_name>
pip install .
```

where `<env_name>` is the name of the environment that will be created.
This only installs dependencies required by the `acdc` packae. You will
need to install anything else manually, e.g.:

```
pip install ipython multiprocessing
```

#### Installing into an existing environment
Activate your desired environment, then `cd` into the cloned repository and install using `pip`:

```
pip install . 
```

## The SQLite Database, `cos_dark`
Currently, two databases are maintained in order to characterize the 
FUV dark rate. One is a SQLite database called `cos_dark.db`- this database
tracks the dark counts over binned regions of the FUV detectors. 
The SQLite database is available on Box- contact Jo Taylor for details.

To access the database, you will need [SQLite](https://www.sqlite.org/index.html) 
if you do not already have it. Once you have it, simply type:

```
sqlite3 cos_dark.db
```

and start exploring! For tips on SQLite queries, click [here](https://www.tutorialspoint.com/sqlite/sqlite_select_query.htm).

Interactive plots of dark counts vs. time are available in the 
`plot_dark_vs_time.ipynb` jupyter notebook.

### SQLite Database Format

The format of the `cos_darks` SQLite database is as follows:

**`darks` table**
| column       | type    | description                                            |
| ------------ | ------- | ------------------------------------------------------ |
| id           | Integer | Primary key ID number                                  |
| filename     | String  | E.g. ipppssoot_corrtag_a.fits                          |
| segment      | String  | COS FUV segment                                        |
| expstart     | Float   | MJD start time                                         |
| exptime      | Float   | Exposure time (s)                                      |
| hv           | Integer | High Voltage                                           |
| latitude     | Float   | Observatory latitude                                   |
| longitude    | Float   | Observatory longitude                                  |
| region       | String  | Region of interest on segment                          |
| region_area  | String  | Area of region                                         |
| xcorr_min    | Integer | Region starting XCORR                                  |
| xcorr_max    | Integer | Region ending XCORR                                    |
| ycorr_min    | Integer | Region staring YCORR                                   |
| ycorr_max    | Integer | Region ending YCORR                                    |
| solar_flux   | Float   | Interpolated solar flux at expstart                    |
| fileloc      | String  | STScI disk location of file                            |
| dark_pha[x]  | Integer | Total dark counts over region for PHA [x] 0 through 31)|

**`solar` table**
| column | type  | description                        |
| ------ | ----- | ---------------------------------- |
| time   | Float | MJD date of solar flux measurement |
| flux   | Float | 10.7cm radio flux                  |

Solar flux values are obtained from [NOAA](https://www.swpc.noaa.gov/phenomena/f107-cm-radio-emissions).

## The MySQL database, `hstcal`
Currently, two databases are maintained in order to characterize the 
FUV dark rate. One is a MySQL database called `hstcal`, hosted on an
STScI internal server.

To access the `hstcal` database, you need MySQL version >= 8. 
For connection information, contact Jo Taylor.

### MySQL Database Format
The format of the `hstcal` MySQL database is as follows. 
There is a table for each HV setting, which includes:
* 163
* 167
* 169
* 171
* 173
* 178
* other

**`darkeventshv<hvsetting>`**

| column       | type        | description                                            |
| ------------ | ----------- | ------------------------------------------------------ |
| id           | Integer     | Primary key ID number                                  |
| xcorr        | Float       | XCORR location of dark event                           |
| ycorr        | Float       | YCORR location of dark event                           |
| pha          | int         | PHA value of dark event                                |
| mjd          | Float(11,5) | MJD time of dark event                                 |
| hv           | Integer     | High Voltage of corresponding dark exposure            |
| segment      | String(5)   | Segment of corresponding dark exposure                 |
| filename     | String(30)  | Filename of corresponding dark exposure                |
| proposid     | Integer     | Proposal ID of corresponding dark exposure             |
