# Copyright (C) 2024 Itay Weintraub ----
# This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
# you are welcome to redistribute it under certain conditions. for details,
# see the GNU General Public License Agreement (in the file COPYING.txt).
rm(list = ls())
# ---------------------------- import statements --------------------------------
suppressPackageStartupMessages({
  suppressWarnings({
    start <- Sys.time()
    require(gridExtra)
    require(deSolve) # solving ordinary differential equations (ODEs)
    require(ggpmisc) # adding statistics to plots
    require(Rcpp) # importing C functions
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(dplyr)
    source("./plotting_functions.R") # various functions for plotting final data
    sourceCpp("./rhs_eval.cpp") # compile external C functions
    source("./input.R")
  })
})
# --------------------------------functions ------------------------------------
# put the results of the numerical integration into a tidy table
organize_data <- function(dat, times, pars,Tenv) {
  dat <- dat %>%
    as.data.frame() %>% # convert to data frame (needed for next step)
    as_tibble() %>% # convert to tibble (tidyverse's improved data frame)
    filter(time %in% times) # only keep specified time points
  names(dat)[1] <- "time" # name the first column "time"
  names(dat)[2] <- paste0("n_", 1, "_", 1) # naming convention:
  names(dat)[3] <- paste0("m_", 1, "_", 1) # naming convention:
  
  dat %>%
    # normalize table by collapsing columns into a key-value column pair
    pivot_longer(cols=2:ncol(.), names_to="variable", values_to="v") %>%
    # split "variable" into value type (density or trait), species, and patch
    separate(variable, c("type", "species", "patch"), sep="_") %>%
    # convert species & patch from string ("1","2",...) to integer (1,2,...)
    mutate(species=as.integer(species), patch=as.integer(patch)) %>%
    # split trait and abundance values into two columns
    pivot_wider(names_from="type", values_from="v") %>%
    mutate(Tenv=Tenv[1:nrow(.)])%>%
    return()
}
# ------------------------------- initial conditions -----------------------------------
ninit <- 1 # reserve memory for initial densities
muinit <- Tmin # initial trait means
ic <- c(ninit, muinit) # merge initial conditions into a vector
pars <- list(rho=rho, kappa=kappa,v=v, nmin=nmin, sigma=sigma,Tmin=Tmin,tE=tE,C=C)
# -------------------------- integrate ODEs -----------------------------------
# consider changing relative (rtol) and absolute (atol) solver tolerances
# for quick/accurate solution trade-off.
#See https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf 
# for elaboration on deSolve package. 
at <-1e-8
rt <-1e-8
maxsteps <- 10000
step <- tE/dpNum
fail_time <- 0
original_tE <- tE
ts <-seq(0, tE, by=step)
Tenv <- ts
for (i in 1:length(ts)) {
  Tenv[i] <- Temp(ts[i],Tmin,C,tE)
}
tryCatch({results <-ode(y=ic, times=seq(0, tE, by=step), func=eqs, parms=pars,
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
    tE <<-floor((fail_time-step)/(step)) * (step) #alternative for round_any
    results <-ode(y=ic, times=seq(0, tE, by=step), func=eqs, parms=pars,
                  method = "bdf",atol  = at, rtol = rt, maxsteps = maxsteps) 
  }
  diagnostics(results)
  dat <- results %>%# put during-climate-change solution into tidy tibble:
    organize_data(times=seq(from=0, to=tE, by=step), pars = pars,Tenv) 
}) 
# --------------------------- generate output ----------------------------------
suppressWarnings(write_csv(dat, path=outfile)) # save data to specified file
pn<-plot_landscape(dat %>% filter(patch %in% c(1)))
pm<-plot_traitLag(dat %>% filter(patch %in% c(1)),nmin)
grid.arrange(pn, pm, ncol=2)
print("R total Runtime")
print(Sys.time()-start)
save.image(file = workspace) # save workspace to specified file