# --------- General parameters---------
S <- 4 # number of species
L <- 20 # number of spatial niches (also called here "patches")
vbar <- 6e-5 # nominal genetic variance (Celsius ^2)
dbar <- 1e-7 # nominal dispersal (1e-7 <=> 1 meter per year)
# more precisely, in units of pole to equator distance , which is ~100,000 km (1e7 meter)

#-------Temperature profile parameters---------
Tmax <- 25.0 # initial temperature at equator
Tmin <- -15 # initial temperature at poles
Cmax <- 23 # projected temperature increase at poles
Cmin <- 11 # projected temperature increase at equator

periodic <- TRUE # TRUE for sinusoidal profile, FALSE for single rising step
cycles <- 3 # period number, namely : T=sin(2*pi*cycles*t)

small <-FALSE # TRUE for short preparation time, FALSE for long preparation time
tstart <- if (small) -1e5 else -1e8 #preparation time before CC (climatic change) onset
tE <- 2e7 # time span of CC

dpNum <- 200 # number of time stamps for tabular results
#----------------Ecological Model parameters----------
v <-c(1.06,1.19,1.23,1.02) * vbar
d <- c(1.49,1.97,1.52,1.98) * dbar
rho <-c(6.85,10.11,1.09,5.2)

# seed <- 11093 # random seed for variations over species properties
# set.seed(seed) # set random seed for reproducibility
# v <- runif(S, 1.0*vbar, 1.25*vbar) # genetic variances
# d <- runif(S, 1.0*dbar, 2.0*dbar) # dispersal rates
# rho <- runif(S, 0.1, 11) # growth-tolerance trade off parameter

kappa <- 0.1 # intrinsic mortality parameter
eta <- 1 # competition width in Celsius
nmin <- 1e-5 # below this threshold density, genetic variances are reduced
aw <- 0.1 # (negative) slope of trait-dependence of tolerance width
bw <- 4 # intercept of trait-dependence of tolerance width
#--------I/O files-----------------
runID <-"trial"
file <- runID
outfile <- paste("outputs/",file, sep = "")
workspace <-paste("parameters/",file, sep="")