runID <-"trial"
tE <-1e6 # simulation time span in years
kappa <- 1e-4 # intrinsic mortality rate
v <- 0.5 # genetic variance in Celsius^2
sigma <- 2 # temperature tolerance
rho <- 2e-4 # growth-tolerance trade off parameter
nmin <- 1e-5 # minimal density threshold , extinction/survival boundary
C <- 10 # total temperature shift in Celsius degrees
T0 <- 15 # initial temperature
dpNum <- 200 # number of time stamps for tabular results

file <- runID
outfile <- paste("outputs/",file, sep = "") 
workspace <-paste("parameters/",file, sep="")