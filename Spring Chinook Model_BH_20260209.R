# =========================
# Spring Chinook LCM based on NOAA's ASRP
# https://bitbucket.org/noaalcm/asrp/src/master/lcm/params/params.spring.chinook.R
# =========================

library(dplyr)
library(tidyverse)
library(ggplot2)

# ------------------------------
# 1. Initial Conditions
# ------------------------------
pop_init <- 3000
n_years  <- 10
n_mc     <- 25
set.seed(123)

# ------------------------------
# 2. Life-History Parameters
# ------------------------------

eggs_per_female <- 5400
prop_female     <- 0.5
egg_capacity    <- 8e6

redd_survival   <- 0.99
emergence       <- 0.40
fry_to_parr_fw  <- 0.30

prop_parr <- 0.5
prop_fry  <- 0.5

# Return Propensity
b2 <- 0.005
b3 <- 0.097
b4 <- 0.60
b5 <- 0.80

Hr    <- 0.056
S_up  <- 0.76
S_pre <- 0.70

# ------------------------------
# 2b. Survival Ranges
# ------------------------------

# Bay survival ranges
bay.parr.min <- 0.04
bay.parr.max <- 0.09
bay.fry.min  <- 0.0005
bay.fry.max  <- 0.002

# Ocean survival ranges by age
so.range <- list(
  c(0.55,0.70),
  c(0.65,0.80),
  c(0.70,0.90),
  c(0.80,0.95),
  c(0.85,0.97)
)

# ------------------------------
# 3. Environmental Stochasticity
# ------------------------------

#---------- Floods ------------#
#currently not being used, NOAA's function for egg
#survival as a function of flood recurrence is plugged
#into the model. 

#flood_probs <- data.frame(
#  period = c("Mid","Late"),
#  typical = c(0.71,0.63),
#  major   = c(0.27,0.34),
#  catastrophic = c(0.016,0.027)
#)

#flood_mults <- c(
#  typical = 0.90,
#  major   = 0.70,
#  catastrophic = 0.40
#)

#get_climate_period <- function(year){
#  if(year <= 35) "Mid" else "Late"
#}

#draw_flood <- function(year){
#  period <- get_climate_period(year)
#  probs  <- flood_probs[flood_probs$period==period,2:4]
#  type   <- sample(names(probs),1,prob=as.numeric(probs))
#  flood_mults[type]
#}

#---------- Temperature ------------#

#temp_probs <- data.frame(
#  period=c("Mid","Late"),
#  normal=c(0.6,0.4),
#  warm=c(0.3,0.4),
#  hot=c(0.1,0.2)
#)

#temp_mults <- c(normal=1.0,warm=0.85,hot=0.70)

#draw_temperature <- function(year){
#  period <- get_climate_period(year)
#  probs  <- temp_probs[temp_probs$period==period,2:4]
#  type   <- sample(names(probs),1,prob=as.numeric(probs))
#  temp_mults[type]
#}

###NOAA temp function
temp_func <- function(t){
  ifelse(t < 18,
         1,
         ifelse(t < 24,
                -(1/6) * t + 4,
                0))
}

#Avg basin-wide temps for mid- and late-century
#based on NOAA's thermalscape_temps data 
#   (just took mean temp from all habitat sections)
mid_temp <- 16.99
late_temp <- 18.74

interp_temp <- function(yr, mid_temp, late_temp, max_year=100){
  k <- 0.08   # steepness
  x0 <- max_year/2
  frac <- 1 / (1 + exp(-k * (yr - x0)))
  mid_temp + (late_temp - mid_temp) * frac
}
# ------------------------------
# 4. Various Functions
# ------------------------------

#Density Dependence
BH.func <- function(S,p,c){
  c[c==0] <- 1
  recruits <- (S*p)/(1+(p/c)*S)
  recruits[recruits<1] <- 0
  recruits
}

# Inverse logit function
# helper function called by another function
inv.logit <- function(x){
  exp(x)/(exp(x)+1)
}

# Change rescale egg survival to 0:1 scale
# helper function called by another function
range.rescale <- function(x,min=0.005554976,max=0.1322917){
  (x-min)/diff(range(c(min,max)))
}
# rescales x the range of predictions,
# where x is the inv.logit from a fitted 
# logit(surv) to recurrence interval relationship:
# logit(surv) ~ -1.88084588 - 0.05511075 * peak.winter.flow
# from Skagit River Age0+ Chinook salmon
# egg-outmigrating juveniles data, rescaled to
# 0 to 1 scale, with the minimum surv of 0.005554976 
# at high winter flows, and
# a maximum surv of 0.1322917 at low winter flows

# Egg survival as a function of flow recurrence
egg.flow.dec <- function(){
  temp.flow <- sample(seq(0.05,100,1),1,
                      prob=1/seq(0.05,100,1))
  range.rescale(inv.logit(-1.88084588-0.05511075*temp.flow))
}

# Bay survival function ----
bay.surv.func <- function(min,max){
  runif(1,min,max)
}


# Marine survival function ----
So.func <- function(min,max){
  runif(1,min,max)
}

# ------------------------------
# 5. Monte Carlo Simulation
# ------------------------------

results_mc <- matrix(NA,nrow=n_mc,ncol=n_years)

for(mc in 1:n_mc){
  
  adults <- pop_init
  
  for(yr in 1:n_years){
    
    females <- adults * prop_female
    
    # ---- Define flow/temp impacts ----
    #flood_mult <- draw_flood(yr)
    #temp_mult  <- draw_temperature(yr)
    trend_temp <- interp_temp(yr, mid_temp, late_temp)
    temp_value <- rnorm(1, trend_temp, 0.5)
    temp_mult  <- temp_func(temp_value)
    flow_mult  <- egg.flow.dec()
    
    # ---- Spawners → Eggs ----
    eggs_raw <- BH.func(females,eggs_per_female,egg_capacity)
    
    eggs_post <- eggs_raw * 
      redd_survival *
      #flood_mult * 
      flow_mult
    
    # ---- Eggs → Fry ----
    fry <- eggs_post * emergence
    
    fry_mig  <- fry * prop_fry
    parr_raw <- fry * prop_parr
    
    parr <- parr_raw * fry_to_parr_fw * temp_mult
    
    # ---- Bay Survival ----
    bay_fry_surv  <- bay.surv.func(bay.fry.min,bay.fry.max)
    bay_parr_surv <- bay.surv.func(bay.parr.min,bay.parr.max)
    
    fry_bay  <- fry_mig * bay_fry_surv
    parr_bay <- parr    * bay_parr_surv
    
    # ---- Ocean Survival Draws ----
    so1 <- So.func(so.range[[1]][1],so.range[[1]][2])
    so2 <- So.func(so.range[[2]][1],so.range[[2]][2])
    so3 <- So.func(so.range[[3]][1],so.range[[3]][2])
    so4 <- So.func(so.range[[4]][1],so.range[[4]][2])
    so5 <- So.func(so.range[[5]][1],so.range[[5]][2])
    
    # ---- Age Structure ----
    age1 <- fry_bay + parr_bay
    age2 <- age1 * so1 * (1-b2)
    age3 <- age2 * so2 * (1-b3)
    age4 <- age3 * so3 * (1-b4)
    age5 <- age4 * so4 * (1-b5)
    age6 <- age5 * so5
    
    returns <- (age1*b2 + age2*b3 +
                  age3*b4 + age4*b5 +
                  age6)
    
    adults_eff <- returns *
      (1-Hr) * S_pre * S_up *
      temp_mult
    
    results_mc[mc,yr] <- adults_eff
    adults <- adults_eff
  }
}

# ------------------------------
# 6. Summaries
# ------------------------------

summary_df <- data.frame(
  Year=1:n_years,
  Mean=apply(results_mc,2,mean),
  SD=apply(results_mc,2,sd),
  P5=apply(results_mc,2,quantile,0.05),
  P95=apply(results_mc,2,quantile,0.95)
)
summary_df
# ------------------------------
# 7. Plot
# ------------------------------

long_mc <- as.data.frame(as.table(results_mc))
colnames(long_mc) <- c("Sim","Year","Adults")
long_mc$Year <- as.numeric(long_mc$Year)

ggplot()+
  geom_line(data=long_mc,
            aes(Year,Adults,group=Sim),
            alpha=0.55,color="steelblue")+
  geom_line(data=summary_df,
            aes(Year,Mean),
            color="darkred",linewidth=1.2)+
  theme_minimal(base_size=14)+
  labs(title="Spring Chinook Adult Returns",
       x="Year",y="Adults")
