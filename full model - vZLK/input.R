# --------- General parameters---------
S <- 4 # number of species
L <- 20 # number of spatial niches (also called here "patches")
vbar <- 6e-5 # nominal genetic variance (Celsius ^2)
dbar <- 1e-3 # nominal dispersal (1e-7 <=> 1 meter per year)
# more precisely, in units of pole to equator distance , which is ~100,000 km (1e7 meter)
#----------------Ecological Model parameters----------
v <- rep(vbar, S) # genetic variances
d <- rep(dbar, S) # genetic variances
seed<-3690
set.seed(seed)
rho <- runif(S, 0.1, 11) # growth-tolerance trade off parameter
kappa <- 0.1 # intrinsic mortality parameter
eta <- 1 # competition width in Celsius
nmin <- 1e-5 # below this threshold density, genetic variances are reduced
aw <- 0.1 # (negative) slope of trait-dependence of tolerance width
bw <- 4 # intercept of trait-dependence of tolerance width
#--------I/O files-----------------
runID <-"run"
file <- runID
outfile <- paste("outputs/",file, sep = "")
workspace <-paste("parameters/",file, sep="")