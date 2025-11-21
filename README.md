# Merlino2025_IL27_ZIKV

[![DOI](https://zenodo.org/badge/989098542.svg)](https://doi.org/10.5281/zenodo.17672080)

## Required software

R >= 4.4.2

##  Installation

Clone the repository

```
git clone https://github.com/amsesk/Merlino2025_IL27_ZIKV.git
```

Install R library dependencies

```
# Change directory to repository
cd Merlino2025_IL27_ZIKV

# Enter R interactive session
R

# Update dependencies from lockfile and exit
> renv::restore(lockfile="renv.lock")
> q()
```

### Generate the manuscript figures

The code for generating the manuscript figures is all contained in `figures.R`. Output folders will be created and populated by the script. 

```
Rscript figures.R
```

**Note**: If you run it from the repository, then no edits to the script are required. I you are running from a different location, then you'll need to update the value of the `REPO` variable in `figures.R` to reflect the location of this repository. 


