# Astro-Ecology-Paper
This rep bears the source code required to reconstruct results shown on "Adaptive Habitability of Exoplanets:Thriving Under Extreme Environmental Change" paper.
The model implemented in rhs_eval.cpp is a slightly thinner version of the code and ecological model used in Akesson et al 2021, https://doi.org/10.1038/s41467-021-24977-x.
Model parameters were stretched to fit geological timescale climatic change, and incorporate other temperature trends (Ideal sinusoidal, vZLK induced).

File Structure:
1.Kozai.R - solves the full octuple order of the vZLK secular dynamics depicted by Naoz 2016. From eccentricity computes temperatures for equator and pole environments. The orbital parameters and resulting temperature profile can be saved and retrieved later.

2.

The computational speed relies on cpp doing the hard work, while the R main script is only there for data preparation and visualization using the ggplot2 plotting functions.
Pay attention that stiff ODE solver updates step size on the go,and so calls for rhs_eval are very frequent for most. 
For that reason, the smaller memory passed to rhs_eval as arguments, the better. 
For pre-computed T configurations like the vZLK, keep the T vector as small as possible. 
For instance , 200 T data points will result in 1-10 minutes run, in accordance to Akesson original performance.  