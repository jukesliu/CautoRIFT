# CautoRIFT
[Jukes Liu](https://github.com/jukesliu). Department of Geosciences and Cryosphere Remote Sensing and Geophysics (CryoGARS) Lab, Boise State University.
#### Contact: jukesliu@boisestate.edu, jukes.liu@gmail.com

<img src="https://github.com/julialiu18/SK-surge-mapping/assets/48999537/8a1748c1-573f-4a30-9618-e7e95c424004" width="600">

200-meter resolution velocity map of Sít' Kusá (Turner Glacier), southeast Alaska, during its 2020 surge produced from a Sentinel-2 image pair from May 5 and May 10 2020. Vectors indicate flow direction.

### Summary
This repository contains code to map glacier surface velocities using the NASA JPL open-source [geogrid and autoRIFT](https://github.com/nasa-jpl/autoRIFT) software using custom feature-tracking parameters. These custom parameters include a reference velocity map, a reference DEM, search distances, chip sizes, date separations, and a stable surface mask (binary mask over all non-moving surfaces). The custom autoRIFT script (`CautoRIFT.ipynb`) may be applied to optical satellite image pairs from PlanetScope, Sentinel-2, and Landsat 7-9. __NOTE: All input images must be cropped to the same extent using an Area of Interest (AOI) shapefile!__ You can manually compile your input satellite images prior to CautoRIFT or the `LS7-9_download_AWS.ipynb` and `S2_COG_download_AWS.ipynb` scripts provided may be used to automatically download and crop the Landsat and Sentinel-2 images, respectively. A separate GitHub repository ([planet_tile2img](https://github.com/CryoGARS-Glaciology/planet_tile2img)) can be used to download and crop the commercial PlanetScope images.

The `generate_stable_surface_mask.ipynb` script can be used to automatically generate the stable surface mask as the inverse of the glacier outline cropped to the Area of Interest (AOI) extent if the glacier outline shapefile and the AOI shapefile are input. The calculation of stable surface error (SSE) `calculate_vmap_SSE.ipynb` uses the stable surface mask to calculate the velocity error associated with each of the velocity maps produced. Additional details are described in the manuscript associated with the citation below. Please use the following citation to attribute the research use of CautoRIFT:

#### Citation
```
Liu, J., Enderlin, E., Bartholomaus, T., Terleth, Y., Mikesell, T., & Beaud, F. (2024). Propagating speedups during quiescence escalate to the 2020–2021 surge of Sít’ Kusá, southeast Alaska. Journal of Glaciology, 1-12. https://doi.org/10.1017/jog.2023.99
```

## Running the code using containers
### Boise State University Borah Users:

CautoRIFT is on a container on Borah. Move your input files onto Borah and activate the container to run the notebooks using

```
module load apptainer/1.2.5
apptainer run /cm/shared/containers/micromamba-jukes.sif jupyter notebook
```

### Docker users:
```
docker pull ghcr.io/bsurc/jukes-micromamba:latest
```

## Installing the environments with micromamba

Two separate python environments must be installed, one for running the CautoRIFT code and one for all other scripts, which we call __preautorift__. You must be particularly careful not to install new packages into the __cautorift__ environment, otherwise you may encounter compatibility issues.

(0)	Install micromamba
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

### Option A: Using the environment files

```
micromamba env create --file cautorift.yaml

micromamba env create --file preautorift.yaml
```

### Option B: Create the environments manually

#### cautorift

(1) Create a new environment named "cautorift" with python 3.8.16
```
micromamba create -n cautorift python=3.8.16
```

(2) Activate the environment. Then, download the following packages into the environment in the following order:
```
micromamba activate cautorift

micromamba install autorift=1.1.0 gdal=3.0.2 opencv=4.5.0 rasterio=1.2.10 notebook matplotlib pandas -c conda-forge

```
(3) Find the correct __autoRIFT.py__ script (within `micromamba/envs/cautorift/`) using the search bar in Finder and find & replace all the instances of __np.bool__ to __bool__ and __np.int__ to __int__. See my __autoRIFT.py__ path below for reference:

```
~/micromamba/envs/cautorift/lib/python3.8/site-packages/autoRIFT/autoRIFT.py
```

(4) De-activate the environment. Only activate and use it for running the `CautoRIFT.ipynb` script.

```
micromamba deactivate
```

#### preautorift

This environment is to be used for running all other scripts in this repository.

(1) Create a new environment named "preautorift" with python 3.9.6
```
micromamba create -n preautorift python=3.9.6
```

(2) Activate the environment. Then, download the following packages into the environment in the following order:
```
micromamba activate preautorift

micromamba install notebook geopandas=0.10.2 gdal=3.3.1 rasterio=1.2.8 boto3=1.34.87 -c conda-forge

```

## Order of operations

Run the scripts in the following order:

(0) `generate_stable_surface_mask.ipynb`

(1) `LS7-9_download_AWS.ipynb` and `S2_COG_download_AWS.ipynb`

At this point you should have separate folders with the Landsat images and/or the Sentinel-2 images and/or PlanetScope images from [planet_tile2img](https://github.com/CryoGARS-Glaciology/planet_tile2img) that are all cropped to the same extent. The cloudy images must be deleted. You should also have a separate folder of desired custom inputs to autoRIFT (reference x-velocities, reference y-velocities, stable surface mask, reference DEMs) all in geotiff format. You will need to provide the paths to these folders in the `CautoRIFT.ipynb` script.

(2) `CautoRIFT.ipynb`

(3) `calculate_vmap_SSE.ipynb`

## COMING SOON: VIDEO DESCRIBING HOW TO SET UP YOUR INPUTS AND RUN THE EACH SCRIPT IN THE REPOSITORY.

## Funding and Acknowledgements

This research is funded by NASA FINESST Award 80NSSC21K1640 and NSF Award ANS1954006.
