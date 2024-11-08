# ismn-drydowns

This repository contains code for analyzing soil moisture drydowns from the [International
Soil Moisture Network (ISMN)](https://ismn.earth/en/) as detailed in the corresponding
manuscript:

> Morgan, B.E., Araki, R., Trugman, A.T., Caylor, K.K. (*in review*). Ecological and hydroclimatic determinants of vegetation water-use strategies. *Nature Ecology & Evolution*.


## Getting started

1. Clone the repository
```bash
$ git clone git@github.com:ecohydro/ismn-drydowns.git
```

2. Create a virtual environment
```bash
$ cd ismn-drydowns
$ conda env create -f environment.yml
$ conda activate ismn
```

In addition to normal python packages, this will install the `drydowns` package,
a custom package for this project. This package is not currently available on PyPI,
but may be in the future. For now, it is installed from the 
[GitHub repository](https://github.com/ecohydro/drydowns).

3. Unzip results file: `data/ismn_results_star_th0_04may.zip`.
```bash
$ cd data
$ unzip ismn_results_star_th0_04may.zip
```
The results file is >100MB and cannot be included in the repository due to size constraints.

4. Download the ISMN data from [International
Soil Moisture Network (ISMN)](https://ismn.earth/en/).

5. Download ancillary data (CHIRPS, dPET, LAI, etc.) from the appropriate sources.

Note: Steps 4 and 5 are only necessary if re-running the entire analysis (`run_ismn.py`).
To plot the figures in the manuscript, the extracted ISMN drydowns and ancillary data located
in the `data` directory are sufficient.

6. Update `config.ini` with the path to the ISMN data.
Note that the config assumes a certain structure for the soil moisture data and ancillary files,
but this can be modified in the `config.ini` file.


## Contents

### Code
There are three types of files in the `code` directory:

1. **Processing scripts + related files**

   - `calc_map.py`: Calculate mean annual precipitation from CHIRPS data.
   - `calc_mean_pet.py`: Calculate mean annual PET from dPET data.
   - `extract_chirps.py`: Extract daily rainfall data from CHIRPS for coordinates in a CSV file (ISMN stations).
   - `extract_pet.py`: Extract daily PET data from dPET for coordinates in a CSV file (ISMN stations).
   - `ismn_ancillary.py`: Create (+ update) ancillary data for ISMN stations (used in analysis), `ismn_ancillary.csv`.
   - `rainfall.py`: Utility functions for calculating rainfall statistics.
   - `separate_appeears.py`: Separate data downloaded from APPEEARS into separate CSV files by ISMN station (if necessary) + recombine. These files are used as ancillary inputs (GPP, LAI, etc.) in the `drydowns` code and/or in `ismn_ancillary.py`.
   - `separate_pet.py`: Separate dPET data into CSV files by ISMN station. These files are used as ancillary inputs in the `drydowns` code.
   - `soil/`: Module for handling ISMN soil moisture data. `soil/station.py` is a wrapper class for the ISMN Python package to make handling the data easier.
   - `utils.py`: Utility functions for processing ISMN data.

2. **`run_ismn.py`**

   This script runs the drydown extraction and requires the [`drydowns`](https://github.com/ecohydro/drydowns) package.  
   The output of this script is the extracted drydowns, which are saved to a CSV file, `ismn_results_star_th0_04may.csv`.

3. **Plotting scripts + related files**

   - `figs.py`: Code for importing the results + plotting Fig. 3 and supplementary figures. Must be run before `ternary.py` and `surfaces.py`.
   - `conceptual.py`: Code for plotting the conceptual figure in the manuscript (Fig. 1).
   - `ternary.py`: Code for plotting the maps in the manuscript (Fig. 2).
   - `surfaces.py`: Code for plotting the surfaces in the manuscript (Figs. 4-5).
  

### Data
The "data" contained in this repository include the extracted ISMN drydowns 
([ismn_results_star_th0_04may.csv](data/ismn_results_star_th0_04may.csv)) and 
ancillary data used in the analysis.



## Contact
Bryn Morgan, [brynmorgan@ucsb.edu](mailto:brynmorgan@ucsb.edu)