# Multi-sensor optimal data merging
This repository documents the method of multi-sensor optimal data merging (MODM)
for optimally combined large number of different observations. 

#USAGE

If you already have two datasets, you just modify the directories in merge_sic.py
and run 

    python merge_sic.py date(e.g. 20220101)

If you have more than two datasets, you can recursively use the merge in merge_sic.py

The method is described at (also as attached as modm_guide.pdf here)

Wang, K. et al. (2020): Multi-sensor data merging of sea ice concentration and 
   thickness. Advances in Polar Science, 31(1): 1-13, doi:10.13679/J.ADVPS.2019.0016.

Also included is a get_amsr2.py, which downloads and interpolate the data to the grid
you desire. The calculations of the uncertainty follows

Spreen et al. (2008): Sea ice remote sensing using AMSR-E 89-GHz channels. J Geophys 
Res-Oceans, 113(C2): C02S03, doi: 10.1029/2005jc003384.


For citing this code, see [![DOI](https://zenodo.org/badge/473923694.svg)](https://zenodo.org/badge/latestdoi/473923694)

