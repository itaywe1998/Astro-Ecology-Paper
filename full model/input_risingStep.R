# --------- General parameters---------
S <- 4 # number of species
L <- 20 # number of spatial niches (also called here "patches")
vbar <- 6e-5 # nominal genetic variance (Celsius ^2)
dbar <- 1e-7 # nominal dispersal (1e-7 <=> 1 meter per year)
# more precisely, in units of pole to equator distance , which is ~100,000 km (1e7 meter)

#-------Temperature profile parameters---------
Tmax <- 15.0 # initial temperature at equator
Tmin <- -25 # initial temperature at poles
Cmax <- 30 # projected temperature increase at poles
Cmin <- 20 # projected temperature increase at equator

periodic <- FALSE # TRUE for sinusoidal profile, FALSE for single rising step
cycles <- -1 # period number, namely : T=sin(2*pi*cycles*t)

small <-FALSE # TRUE for short preparation time, FALSE for long preparation time
tstart <- if (small) -1e5 else -1e8 #preparation time before CC (climatic change) onset
tE <- 2e6 # time span of CC

dpNum <- 200 # number of time stamps for tabular results
#----------------Ecological Model parameters----------
v <-c(1.1,1.05,1.075,1.075) * vbar
d <- c(1.03,1.01,1.07,1.12) * dbar
rho <-c(6.41,9.05,4.05,2.17)

# seed <- 7234 # random seed for variations over species properties
# set.seed(seed) # set random seed for reproducibility
# v <- runif(S, 1.0*vbar, 2.0*vbar) # genetic variances
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