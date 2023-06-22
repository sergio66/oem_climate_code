## see sergio_matlab_jobB.sbatch : this is to make the rtp files

## 5 = ERA5, 6 = CMIP6, 7 = AIRS L3, 8 = <MERRA2, 9 = AMIP6, 10 = CLIMCAPS, 11 = UMBC
sbatch  --array=1-64 sergio_matlab_jobB.sbatch 5
sbatch  --array=1-64 sergio_matlab_jobB.sbatch 6 
sbatch  --array=1-64 sergio_matlab_jobB.sbatch 7
sbatch  --array=1-64 sergio_matlab_jobB.sbatch 8
sbatch  --array=1-64 sergio_matlab_jobB.sbatch 9 
sbatch  --array=1-64 sergio_matlab_jobB.sbatch 10
sbatch  --array=1-64 sergio_matlab_jobB.sbatch 11

## this is to make the spectral trends, based on whether you run 5,6,7,8,9,10
## sbatch  --array=1-64 sergio_matlab_jobB.sbatch 12
