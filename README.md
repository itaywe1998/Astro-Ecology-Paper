# Astro-Ecology-Paper

This rep bears the source code required to reconstruct results shown on
"Adaptive Habitability of Exoplanets:Thriving Under Extreme
Environmental Change" paper.

The model implemented is a slightly thinner version of the code and
ecological model used in Akesson et al 2021,
<https://doi.org/10.1038/s41467-021-24977-x>. Model parameters were
stretched to fit geological timescale climatic change, and incorporate
other temperature trends (Ideal sinusoidal, vZLK induced).

Repository File Structure: [d-directory, f-file]

1.  (d) simple model

    Simple model described in paper. Assumes single species in a single
    habitat only intrinsic growth and basic directional selection.

2.  (d) full model - Ideal

    Full ecological model,including multiple species and
    habitats,migration and competition. Code adapted to ideal climatic
    changes implementations , namely step function and sinusoidal.

3.  (d) full model - vZLK

    Full ecological model, yet this time opted to accept pre-computed
    temperature table,namely the vZLK secular dynamic induced climatic
    change.

4.  (f) COPYING.txt

    GNU GENERAL PUBLIC LICENSE.

5.  (f) README.md

    this.

Each of the directories 1,2,3 hold specific versions of the same file
structure:

1.  (f) main.R

    Imports inputs, supporting packages, model source code and plotting
    utilities. Runs model using an ODE solver (deSolve) over given time
    span. Summarize and plots results.

2.   (f) input.R

    Incorporates ecological model parameters, I/O file names,
    temperature profile parameters (except for vZKL).

3.  (f) model.cpp

    Contains main step to be run by the ODE solver, called "eqs", in
    addition to temperature profile utilities. Critical for total run
    time, and so implemented in cpp.

4.   (f) plotting_functions.R

    Plotting utilities, e.g. for population density over continuous time
    (plot_landscape) or continuous space (plot_timeseries).

5.  (d) outputs

    Directory for automatic documentation of tabular results over
    defined time series.

6.   (d) parameters

    Directory or automatic documentation of workspace in use.

Exceptions:

-   vZLK using the same plotting file as full model-ideal

-   vZLK in addition contains:

    -   (f) Kozai.R

        Solves the full octuple order of the vZLK secular dynamics
        depicted by Naoz 2016. From eccentricity computes the
        temperatures for equator and pole environments.

    -   (d) kozai_parameters

        Directory to store astrophysical parameters (masses, orbit)
        used.

-   Full model - ideal has 2 example input files:

    -   Input.R - configuration used in paper for sinusoidal climatic
        change

    -   Input_risingStep.R - configuration used in paper for step
        climatic change

**Note:** The computational speed relies on cpp doing the hard work,
while the R main script is only there for data preparation and
visualization using the ggplot2 plotting functions.

Pay attention that stiff ODE solver updates step size on the go,and so
calls for the eqs function in model.cpp are very frequent for most. For
that reason, the smaller memory passed to eqs as arguments, the better.

For pre-computed T configurations like the vZLK, keep the T vector as
small as possible. For instance , 200 T data points will result in 1-10
minutes run, in resemblance to Akesson's original performance.
