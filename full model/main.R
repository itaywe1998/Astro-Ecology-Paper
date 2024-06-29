# Copyright (C) 2024 Itay Weintraub ----
# This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
# you are welcome to redistribute it under certain conditions. for details,
# see the GNU General Public License Agreement (in the file COPYING.txt).
rm(list = ls())
setwd("~/Astro-Ecology-Paper/full model")
# ---------------------------- import statements --------------------------------
suppressPackageStartupMessages({
  suppressWarnings({
    start <- Sys.time()
    require(gridExtra)
    require(deSolve) # solving ordinary differential equations (ODEs)
    require(tidyverse) # manipulating and visualizing data
    require(ggpmisc) # adding statistics to plots
    require(Rcpp) # importing C functions
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(dplyr)
    sourceCpp("model.cpp") # compile external C functions
    source("input.R")
    source("plotting_functions.R") # various functions for plotting final data
  })
})
# --------------------------------functions ------------------------------------
# put the results of the numerical integration into a tidy table
organize_data <- function(dat, times, pars) {
  dat <- dat %>%
    as.data.frame() %>% # convert to data frame (needed for next step)
    as_tibble() %>% # convert to tibble (tidyverse's improved data frame)
    filter(time %in% times) # only keep specified time points
  names(dat)[1] <- "time" # name the first column "time"
  index <- 2 # keep track of which column we are naming with this counter
  for (k in 1:pars$L) {
    for (i in 1:pars$S) { # name columns for densities
      names(dat)[index] <- paste0("n_", i, "_", k) # naming convention:
      index <- index + 1 # "type_species_patch" - type is either m (trait),
    } # or n (density)
  }
  for (k in 1:pars$L) {
    for (i in 1:pars$S) { # name columns for trait values
      names(dat)[index] <- paste0("m_", i, "_", k) # (same naming convention)
      index <- index + 1
    }
  }
  dat %>%
    # normalize table by collapsing columns into a key-value column pair
    pivot_longer(cols=2:ncol(.), names_to="variable", values_to="v") %>%
    # split "variable" into value type (density or trait), species, and patch
    separate(variable, c("type", "species", "patch"), sep="_") %>%
    # convert species & patch from string ("1","2",...) to integer (1,2,...)
    mutate(species=as.integer(species), patch=as.integer(patch)) %>%
    # split trait and abundance values into two columns
    pivot_wider(names_from="type", values_from="v") %>%
    return()
}
# ------------------------------- initial conditions -----------------------------------
# dispersal matrix
mig <- matrix(0, L, L) # initialize dispersal matrix
for (k in 2:L) mig[k-1,k] <- 1 # each species can only migrate to the two
mig <- mig + t(mig) # nearest-neighbor patches
ninit <- matrix(0, S, L) # reserve memory for initial densities
muinit <- matrix(seq(Tmin, Tmin, l=S), S, L) # initial trait means
Tempinit <- Temp(seq(from=0, to=1, l=L), 0, tE, Cmax, Cmin, Tmax, Tmin, periodic, cycles)
for (i in 1:S) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
ic <- c(ninit, muinit) # merge initial conditions into a vector
# coerce parameters into a list
pars <- list(S=S, L=L, rho=rho, kappa=kappa, eta=eta,
             vbar=vbar,v=v, nmin=nmin, aw=aw, bw=bw,
             Tmax=Tmax, Tmin=Tmin, Cmax=Cmax, Cmin=Cmin,
             tE=tE, d=d, mig=mig, periodic=periodic, cycles=cycles)
# --------------------------- integrate ODEs -----------------------------------
# consider changing relative (rtol) and absolute (atol) solver tolerances
# for quick/accurate solution trade-off.
#See https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf 
# for elaboration on deSolve package. 
at <-1e-10
rt <-1e-10
before_step <- -tstart/(dpNum)
tryCatch({before_cc <-ode(y=ic, times=seq(tstart, 0, by=before_step), func=eqs, parms=pars,
                          method="bdf", atol  = at, rtol = rt, maxsteps = 10000)},
         error=function(e){message("All Species Extinct")
           return(NA)}) # integrate ODEs before climate change starts
diagnostics(before_cc)
ic <- as.numeric(before_cc[nrow(before_cc),-1]) # final state -> new initial cond.
before_cc <- before_cc %>% # put before-climate-change solution into tidy tibble:
  organize_data(times=seq(from=tstart, to=0, by=before_step), pars = pars) %>%
  filter(time!=0) # remove time point 0 (will be starting point of during_cc)
print("Before CC")
during_step <- tE/dpNum
at <-1e-10
rt <-1e-10
fail_time <- 0
original_tE <- tE
tryCatch({during_cc <-ode(y=ic, times=seq(0, tE, by=during_step), func=eqs, parms=pars,
                          method = "bdf",atol  = at, rtol = rt, maxsteps = 10000)},
         error=function(fail_time){
           message("All Species Extinct")
           fail_time<<-as.numeric(fail_time$message)},
         finally = {
           if (fail_time > 0) {
             outfile <<- paste(outfile,"_FAILED",sep="")
             unlink(workspace) # Deleting old name workspace
             workspace <<- paste(workspace,"_FAILED",sep="")
             save.image(file = workspace)
             during_step <<- 1000
             tE <<-floor((fail_time-during_step)/during_step) * during_step #alternative for round_any
             # if needed in another place will move to a function
             during_cc <-ode(y=ic, times=seq(0, tE, by=during_step), func=eqs, parms=pars,
                             method = "bdf",atol  = at, rtol = rt, maxsteps = 10000) 
           }
           diagnostics(during_cc)
           during_cc <- during_cc %>% # put during-climate-change solution into tidy tibble:
             organize_data(times=seq(from=0, to=tE, by=during_step), pars = pars) #%>%
           
           # merge data from before, during, and after climate change
           dat <- bind_rows(before_cc, during_cc) %>%
             # merge average genetic variance and dispersal into a single column
             mutate(parameterization=paste0("V=", vbar, " d=", dbar)) %>%
             # create regions
             mutate(region=case_when(
               (patch<=round(max(patch)/3))   ~ "polar", # top third of patches are "polar"
               (patch>=round(2*max(patch)/3)) ~ "tropical", # bottom third are "tropical"
               TRUE                           ~ "temperate")) # the rest are "temperate"
           
         })  # integrate from start to end of climate change
# --------------------------- generate output ----------------------------------
suppressWarnings(write_csv(dat, path=outfile)) # save data to specified file
#plot_timeseries(dat %>% filter(time %in% c(tstart,seq(0, 2e6, by=5e5))))
plot_timeseries(dat %>% filter(time %in% c(tstart,seq(0, 2e7, by=2e6))))
#plot_landscape(during_cc %>% filter(patch %in% c(1,11,20)))
print("Final Runtime")
print(Sys.time()-start)
save.image(file = workspace)