# Copyright (C) 2024 Itay Weintraub ----
# This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
# you are welcome to redistribute it under certain conditions. for details,
# see the GNU General Public License Agreement (in the file COPYING.txt).
rm(list = ls())
setwd("~/Astro-Ecology-Paper/full model - vZLK")
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
    source("kozai.R")
    source("~/Astro-Ecology-Paper/full model - Ideal/plotting_functions.R") # various functions for plotting final data
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
# ---- Kozai ------
# envoke the kozai temprature calculation, unless the workspace already exist
old_profile <- TRUE
if (old_profile){
  wksp_name <- "trail_light"
  kozai_wksp <- paste("kozai parameters/",wksp_name, sep="")
  tmp.env <- new.env() # create a temporary environment
  load(kozai_wksp, envir=tmp.env) # load workspace into temporary environment
  T_kozai <- tmp.env$Tvec
  step <- tmp.env$step / tmp.env$yr
  crit_diff <-tmp.env$indicating_diff / step 
  min_times <- tmp.env$min_times
  max_times <-tmp.env$max_times
  rm(tmp.env) 
}else{
  T_kozai <- kozai()
}
# ------------------------------- initial conditions -----------------------------------
lT <- length(T_kozai[,1])
Tmin <- unname(T_kozai[1,2])
Tmax <- unname(T_kozai[1,3])
# dispersal matrix
mig <- matrix(0, L, L) # initialize dispersal matrix
for (k in 2:L) mig[k-1,k] <- 1 # each species can only migrate to the two
mig <- mig + t(mig) # nearest-neighbor patches
ninit <- matrix(0, S, L) # reserve memory for initial densities
muinit <- matrix(seq(Tmin, Tmin, l=S), S, L) # initial trait means
Tempinit <- Temp(seq(from=0, to=1, l=L), Tmax, Tmin)
for (i in 1:S) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
ic <- c(ninit, muinit) # merge initial conditions into a vector
# coerce parameters into a list
pars <- list(S=S, L=L, rho=rho, kappa=kappa, eta=eta,
             vbar=vbar,v=v, nmin=nmin, aw=aw, bw=bw,
             d=d, mig=mig, T_kozai=T_kozai, lT=lT)
# --------------------------- integrate ODEs -----------------------------------
# consider changing relative (rtol) and absolute (atol) solver tolerances
# for quick/accurate solution trade-off.
#See https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf 
# for elaboration on deSolve package. 
at <-1e-8
rt <-1e-8
maxsteps <- 10000
tE <-tail(T_kozai, n=1)[1]
step <- (unname(T_kozai[2,1])-unname(T_kozai[1,1]))
fail_time <- 0
original_tE <- tE


tryCatch({results <-ode(y=ic, times=seq(0, tE-step, by=step), func=eqs, parms=pars,
                        method = "bdf",atol  = at, rtol = rt, maxsteps = maxsteps)
},
error=function(this){
  message(this$message)
  fail_time<<-as.numeric(this$message)
},
finally = {
  if (fail_time > 0) {
    outfile <<- paste(outfile,"_FAILED",sep="")
    unlink(workspace) # Deleting old name workspace
    workspace <<- paste(workspace,"_FAILED",sep="")
    save.image(file = workspace)
    # lets try without the fail_time - step, to see better what happens
    tE <<-floor((fail_time-step)/(step)) * (step) #alternative for round_any
    # if needed in another place will move to a function
    results <-ode(y=ic, times=seq(0, tE, by=step), func=eqs, parms=pars,
                  method = "bdf",atol  = at, rtol = rt, maxsteps = maxsteps) 
  }

  diagnostics(results)
  results <- results %>% # put during-climate-change solution into tidy tibble:
    organize_data(times=seq(from=0, to=tE, by=step), pars = pars) #%>%
  
  temperature<-rep(seq(Tmin, Tmax, l=L),each = S)
  limit <- nrow(results)/(S*L)
  for (i in 2:limit){
    Tlow <- unname(T_kozai[i,2])
    Thigh <- unname(T_kozai[i,3])
    temperature<-c(temperature, rep(seq(Tlow, Thigh, l=L),each = S))
  }
  dat <-results%>%mutate(Tenv=temperature)
}) 
# --------------------------- generate output ----------------------------------
suppressWarnings(write_csv(dat, path=outfile)) # save data to specified file
times <-c(min_times[c(TRUE,FALSE)],max_times[c(TRUE,FALSE)])
plot_timeseries(dat %>% filter(time %in% c(0,times[c(TRUE,FALSE)])))
#plot_landscape(dat %>% filter(patch %in% c(1,11,20)))
print("Final Runtime")
print(Sys.time()-start)
save.image(file = workspace)