suppressPackageStartupMessages({
  suppressWarnings({
    rm(list = ls())
    start <- Sys.time()
    require(deSolve) # solving ordinary differential equations (ODEs)
    library(ggplot2)
    library(patchwork)
    library(pracma)
  })
})
# --- Functions ----
present_profile <- function(name){
  folder <-"/home/itay/EcoEvolution/Kozai_parameters/"
  file <- paste(folder,name,sep = "")
  dat <- load(file)
  p <-p1+p2+p3+plot_layout(ncol=1)
  print(p)
  save <- FALSE
  if(save){
    plot_path <- paste(folder,"plots/",name,sep = "")
    ggsave(filename =  paste(plot_path,".png",sep ="")
           , plot = p, width = 6, height = 4, units = "in")
  }

}

# Full octuple order of the vZLK secular dynamic. See supplementary materials for Naoz 2016.
kozai_osc <- function(t, state, params){
  # global parameters
  L1 <- unlist(params[1])
  L2 <- unlist(params[2])
  Gtot <- unlist(params[3])
  C3_noG <- unlist(params[4])
  C2_coeff <-unlist(params[5])
  
  # state variables
  e1 <- state[1]
  e2 <- state[2]
  G1 <- state[3]
  G2 <- state[4]
  cosi1 <- state[5]
  cosi2 <- state[6]
  omega1 <- state[7]
  omega2 <-state[8]
  
  # useful parameters per iteration
  C3 <- C3_noG / G2^5
  C2 <- (e2 ^2 -1) * C3 / C2_coeff
  i1 <- acos(cosi1)
  i2 <- acos(cosi2)
  itot <- i1 + i2
  
  cit <- cos(itot)
  co1 <- cos(omega1)
  co2 <- cos(omega2)
  c2o1 <- cos(2*omega1)
  s2o1 <- sin(2*omega1)
  e12 <- e1^2
  e22 <-e2^2
  so1 <- sin(omega1)
  so2 <- sin(omega2)
  sit <- sin(itot)
  
  B <- 2 +5*e12 - 7 * e12 * c2o1
  A <- 4 + 3*e12 -2.5 * B * sit^2
  cosphi <- -co1 * co2 - cit*so1*so2
  
  # Rates
  
  omega1_dt <-6* C2*(1/G1 * (4*cit^2 + (5*c2o1-1)*(1-e12-cit^2)+cit/G2 * (2+e12* (3-5*c2o1)))) -
              C3*e2*(e1*(1/G2 + cit/G1)* (so1 * so2 * (10 * (3*cit^2 - 1)*(1-e12)+A) - 5*B * cit * cosphi)
                       -(1-e12)/(e1 * G1) * (so1*so2 *10 * cit * sit^2 * (1-3*e12) + cosphi*(3*A-10*cit^2 +2)))
  
  omega2_dt <- 3*C2*(2*cit/G1 * (2+e12*(3-5*c2o1)) + 1/G2 * (4+6*e12+(5*cit^2-3)*(2+e12*(3-5*c2o1)))) +
               C3*e1*(so1*so2*((4*e22+1)/(e2*G2) * 10 * cit * sit^2 *(1-e12)  -
                      e2*(1/G1 + cit/G2)*(A+10*(3*cit^2-1)*(1-e12))) +
                      cosphi*(5*B*cit*e2*(1/G1 +cit/G2) + (4*e22+1)/(e2*G2) * A)  
                        )
  
  e1_dt <- C2 * (1-e12)/G1 * (30*e1*sit^2 * s2o1) + C3*e2*(1-e12)/G1 * (35*cosphi*sit^2 *e12 * s2o1
            -10*cit*sit^2 * co1 * so2 * (1-e12) - A*(so1*co2-cit*co1*so2))
  e2_dt <- -C3*e1*(1-e12)/G2 * (10*cit*sit^2 * (1-e12)*so1*co2+ A*(co1*so2-cit*so1*co2))
  G1_dt <- -C2*30*e12*s2o1*sit^2 +C3*e1*e2*(-35*e12*sit^2 *s2o1 * cosphi +A*(so1*co2-cit*co1*so2) + 10*cit*sit^2 * (1-e12)*co1*so2)
  G2_dt <- C3 * e1*e2*(A*(co1*so2-cit*so1*co2) +10*cit*sit^2 * (1-e12)*so1*co2)
  ### Not a part of the state vector but yet useful intermediate equations
  H1_dt <- G1/Gtot * G1_dt - G2/Gtot * G2_dt
  H2_dt <- -H1_dt
  ####
  cosi1_dt <- H1_dt/G1 - G1_dt/G1 * cosi1
  cosi2_dt <- H2_dt/G2 - G2_dt/G2 * cosi2
  
  
  list(c(e1_dt,e2_dt,G1_dt,G2_dt,cosi1_dt,cosi2_dt,omega1_dt,omega2_dt))
}

lyr_to_AU <- function(l){
  return(l*6.3241e4)
}

# Convert mili-arcseconds to distance in the same units as the radius
mas_to_dist <- function(angle, r){
  p<-2*pi*r
  as <- (p/360)/(60*60)
  mas <- 1e-3 * as
  return(angle*mas)
}

to_radians<-function(deg){
  return(deg*pi/180)
}

to_deg<-function(rad){
  return(rad*180/pi)
}



#* A little "Controls & sensitivities" Guide
#*   1.if desired to decrease the gap between the pole and equator - increase obliquity (eps)
#*   even possible to revert the relation above eps=45 degrees, making the pole hotter than the equator
#*   2. the Jupiter-like perturbuter's mass (m3) affects the frequency of the eccentricity (e1)
#*    certain values will generate a periodic (ordered) pattern, others will create a rising frequency
#*   3. i1 (inclination for planet) heavily impacts - increasing makes the eccentricity much more extreme and hence the temperature difference over time
#*   4. lowering i2 as well depresses to a well ordered and low difference eccentricity, so it is the combination of both inclinations who rules  
#*   5. Additionally, i1-i2 matters, the larger it is the less moderate the system (e1) is
#*   6. As planets' mass tend to be lower compared to the other masses, it is more vulnerable to extreme changes
#*   the bigger it is the more ordered behavior occurs. In small masses - rising frequency
#*   7. Omegas don't change much (as expected), just the initial point of the same behavior
#*   8. initial e1 changes the stability of the e1 solution, leading to various resulting 
#*   trends, not very linearly responding.
#*   9. Star mass does not change much of its own, but in order to keep sense,
#*   the mass-luminosity relations obliges to increase much of Ls with minor changes to Ms
#*   10. In a similar manner - a1 also controls the average heat, similarly to Ls
#*   and so both parameters are recommended to stay the same.
#*   11. Albedo is a rather straight forward tuning to average heat 
#*   12. the more distant (a2) the perturbater is the more ordered behavior 
#*   extremely close Jupiter-like breaks the solution (as expected)
#*
#*
#*
#*
#*
#*   
kozai <-function(){
  # -----Global Constants-----
  # working in mks
  AU <- 1.495978707e11 # m
  yr <- 31556926 # s
  Ls <- 3.846e26 # W
  Ms <- 1.989e30 # kg
  Mj <- 1.89813e27 # kg
  Me <- 5.972e24 #kg
  GC <- 6.6743e-11
  sigma <- 5.670373e-8 # Stephan-Boltzmann constant [W m^-2 K^-4]
  albedo <- 0.29 # Earth's Albedo
  eps <- to_radians(24) # Earth's Obliquity
  T0 <- -273.15 # Kelvin to Celsius conversion
  
  #---System data ------
  m1 <- 1.000 * Ms # Star
  m2 <- 1 * Me # Terrestrial planet
  m3 <- 7.5* Mj # Gas-giant , can change back to 1 to see more dense repetitions
  # 1 is for inner binary
  a1 <- 0.51 * AU
  e1 <- 	0.01
  i1 <- to_radians(64.9)
  omega1 <- to_radians(0)
  # 2 is for outer binary
  a2 <- 50* AU
  e2 <- 0.01
  i2 <- to_radians(0.1)
  omega2 <- to_radians(0) 
  
  rel <- a1/a2
  #-------Kozai Inputs-----------
  # constants
  L1 <-m1*m2/(m1+m2) * (GC * (m1+m2) * a1)^0.5
  L2 <-m3*(m1+m2)/(m1+m2+m3) * (GC * (m1+m2+m3) * a2)^0.5
  # Do not change order of terms since in certain combinations, using mks with the relevant masses will 
  # overflow the intermediate result either smaller than 1e-300 or bigger than 1e+300, which are the
  # numerical boundaries in R. In that case, following computations will be corrupted by NaN.
  C3_noG <-  -15/16 * GC^2 /4 * (m1+m2)^9 / (m1+m2+m3)^4 * (m1-m2)/(m1*m2)^5 * L1^6 / L2^3  * m3^9
  C2_coeff <- 15/4 * (m1-m2)/(m1+m2) * (a1/a2)
  # initial values for state vector
  G1 <-L1 * (1-e1^2)^0.5
  G2 <-L2 * (1-e2^2)^0.5
  H1 <- G1 * cos(i1)
  H2 <- G2 * cos(i2) 
  # constant
  Gtot <- H1 + H2
  
  ic <- c(e1, e2, G1, G2, cos(i1), cos(i2), omega1, omega2)
  pars <- list(L1 = L1 , L2 = L2, Gtot = Gtot, C3_noG = C3_noG , C2_coeff = C2_coeff)
  at <- 1e-10
  rt <- 1e-10
  tE <- 1.001e9 * yr
  stepNum <- 1002
  step <- tE/stepNum
  
  workspace_name <- "KozaiPreciseDesignForPaper" #design3 is extra points at the end so failtime=1e9
  workspace <- paste("~/EcoEvolution/Kozai_parameters/betterRes_",workspace_name, sep="")
  
  #---- Differential Equation -------
  results <-ode(y=ic, times=seq(0, tE, by=step), func=kozai_osc, parms=pars,
                method="bdf", atol  = at, rtol = rt, maxsteps = 5000)
  diagnostics(results)
  
  sampled_results <- results[seq(0, stepNum, by=1), ]
  disp_results <- sampled_results
  disp_results[,1] <- sampled_results[,1]/yr
  disp_results[,6:7] <-to_deg(acos(sampled_results[,6:7]))
  disp_results[,8:9] <-sampled_results[,8:9] %% 360
  times <- disp_results[,1]
  ecc_vec <- disp_results[,2] #e1
  disp_results <-as.data.frame(disp_results)
  colnames(disp_results) <- c("Time(years)", "e1", "e2", "G1", "G2", "i1", "i2", "omega1", "omega2")
  
  p1<-ggplot(data = as.data.frame(cbind(times,ecc_vec)))+aes(x= times, y= ecc_vec) + ylab("e")+xlab("time (years)")+
    geom_line()+ theme(axis.text= element_text(size = 22),axis.title = element_text(size=24))
  
  
  
  
  #----- Eccentricity to Temperature ----------
  d <- a1
  Lum <- 1 * Ls 
  e <- ecc_vec # This should come from kozai results
  # sources :
  #basic calculation :https://web.archive.org/web/20210605120431/https://scied.ucar.edu/earth-system/planetary-energy-balance-temperature-calculate
  # eccentricity and obliquity addition: https://www.sciencedirect.com/science/article/pii/S1631071310000052
  # the relation of average energy received over entire orbit to eccentricity
  E <- (Lum*(1-albedo)/(16*sigma*pi*d^2))*(1-e^2)^(-0.5) 
  E_pole <- E * sin(eps)/pi
  E_equator <- E * 2/pi^2 * ellipke(sin(eps))$e
  T_pole <- (E_pole)^(1/4) + T0
  T_equator <- (E_equator)^(1/4) + T0
  
  p2<-ggplot(data = as.data.frame(cbind(times,T_equator)))+aes(x= times, y= T_equator) + ylab(expression('T'[equator]))+xlab("time (years)")+
    geom_line()+ theme(axis.text= element_text(size = 22),axis.title = element_text(size=24))
  
  p3<-ggplot(data = as.data.frame(cbind(times,T_pole)))+aes(x= times, y= T_pole) + ylab(expression('T'[pole]))+xlab("time (years)")+
    geom_line()+ theme(axis.text= element_text(size = 22),axis.title = element_text(size=24))
  
  p1+p2+p3+plot_layout(ncol=1)
  
  Tvec <-cbind(times,T_pole, T_equator)
  delta_equator <- max(T_equator)-min(T_equator)
  delta_pole <- max(T_pole) - min(T_pole)
  
  Tgap <- T_equator-T_pole
  avg_Tgap <- mean(Tgap)
  min_Tgap <- min(Tgap)
  max_Tgap <- max(Tgap)
  
  indicating_diff <- max(abs(diff(T_equator))) # the bottleneck of evolution
  # for a single step - the biggest temperature difference
  ecc_diff <- max(abs(diff(ecc_vec)))
  
  max_indices <- which(diff(sign(diff(ecc_vec)))==-2)+1
  min_indices <- which(diff(sign(diff(ecc_vec)))==2)+1
  max_times <- times[max_indices]
  min_times <- times[min_indices]
  
  save.image(file = workspace)
  
  Tvec
  indicating_diff
  
}
