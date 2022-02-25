# COS Dark Characterization

## Another background correction?
COS spectroscopic science in the extreme UV regime is limited by detector 
background noise. The COS pipeline, `calcos` does not perform an optimized 
background subtraction. In order to achieve the maximum scientific value of
the COS instrument, we have a designed a custom characterization and correction
of the COS FUV dark rate.

This package creates and maintains databases needed to measure the
dark rate as a function of time, position, and PHA. It also creates
FUV superdarks and uses them to perform custom dark corrections to 
COS FUV data.

## Installation
If you do not already have Conda installed, you need to download and install
either Miniconda or Anaconda. Miniconda provides a bare minimum Conda
environment. Anaconda provides a full Conda root environment along with
many other tools, libraries, and utilities.
* get [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
* get [Anaconda](https://www.anaconda.com/products/individual)

Now clone this repository on your local machine. Do this by typing:
`git clone git@github.com:spacetelescope/cos_dark.git`
for SSH connection, or:
`git clone https://github.com/spacetelescope/cos_dark.git`
for HTTPS connection.

### Installing into a fresh environment
Go to the cloned repository, and from a bash shell enter:

```
conda create -n <env_name> python
conda activate <env_name>
pip install -e .
```

This only installs packages needed for the `cos_dark` software. You will
need to install anything else manually, e.g.:

```
pip install ipython scipy asdf
```

### Installing into an existing environment
Go to the cloned directory and from a bash shell enter:

```
pip install -e . 
```

## COS Dark Databases
Currently, two databases are maintained in order to characterize the 
FUV dark rate. One is a SQLite database called `cos_dark.db`- this database
tracks the dark counts over binned regions of the FUV detectors. 
The other is a MySQL database called `hstcal`, hosted on the
STScI internal server `PLTANMYSQL`. 

## The SQLite Database, `cos_dark`

If you are within the STScI network you can create a local version
of the the SQLite database like so:

```
python create_db.py
python update_db.py
```

If you are outside the STScI network, you can use the sqlite3 database 
available on Box- contact Jo Taylor for details.

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

**Darks**
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

**Solar**
| column | type  | description                        |
| ------ | ----- | ---------------------------------- |
| time   | Float | MJD date of solar flux measurement |
| flux   | Float | 10.7cm radio flux                  |

Solar flux values are obtained from [NOAA](https://www.swpc.noaa.gov/phenomena/f107-cm-radio-emissions).

## The MySQL database, `hstcal`

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

**darkeventshv\<hvsetting>**

| column       | type    | description                                            |
| ------------ | ------- | ------------------------------------------------------ |
| id           | Integer | Primary key ID number                                  |
| xcorr        | Float  | XCORR location of dark event                            |
| ycorr        | Float  | YCORR location of dark event                            |
| pha          | int   | PHA value of dark event                                  |
| mjd          | Float(11,5)   | MJD time of dark event                           |
| hv           | Integer | High Voltage of corresponding dark exposure            |
| segment      | String(5)   | Segment of corresponding dark exposure             |
| filename     | String(30)   | Filename of corresponding dark exposure           |
| proposid     | Integer  | Proposal ID of corresponding dark exposure            |
