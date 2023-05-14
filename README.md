# 21cm_beam
Primary beam measurement for 21cm interferometers

This repository contains the scripts using the direct optimal mapping ([Xu et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...938..128X/abstract)) to measure the primary beam.

The strategy is to trace the amplitude of one specific source as it transits. The relative change in amplitude measures the profile of the array-averaged primary beam at one specific slice.
With many sources, we measure the primary beam at different slices to eventually reconstruct the 2D primary beam map.

The direct optimal mapping can map a specific part of the sky regardless where we are pointing. 
We only compute pixels around the source of interest to minimize computation. The point-spread-function matrix (P matrix), more details can be find in [Xu et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...938..128X/abstract), is used to estimate power leakage from other bright sources to our source of interest via grating lobes.
