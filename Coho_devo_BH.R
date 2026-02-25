#==================================
# 
# 
#  Coho LCM based on SEPA LCM utilizes EDT outputs:
#  productivity
#  capacity
#  prespawner survival
#  juvenile to spawner survival
# 
# USING KS EDT inputs
#
#=================================

#Bring in data

dat <- read.csv("C:/Users/Blake.Hamilton/OneDrive - Kleinschmidt Associates/Projects/Chehalis/Pop ass model/EDT Outputs/Master EDT Outputs.csv",
                check.names = FALSE) %>% 
  filter(Species == 'Coho')

#Filter for only subpops of interest
dat <- dat %>%
  filter(
    Subpopulation %in% c("Above Crim SB", "Chehalis RBF to Crim SB")
  )

#===  Functions  

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

g <- function(x) head(x)

#=== Set up loop objects
subpops   <- unique(dat$Subpopulation)
scenarios <- unique(dat$Scenario)

results_eq_coho <- data.frame()


# ================================
# COHO SALMON LCM (Figure 3)
# ================================

beverton_holt <- function(N, p, c) {
  (p * N) / (1 + (p / c) * N)
}

make_A_marine_coho <- function(
    SAR,
    prespawn_surv,
    ocean_surv = c(0.7, 0.7),   # Table 2
    b2 = 0.033                  # age-2 maturation (jacks)
) {
  
  # Partition Grays Harbor survival from SAR (QA/QC section p.24)
  gh_surv <- SAR / prod(ocean_surv)
  
  # State order:
  # 1 Spawners
  # 2 Parr_end_summer
  # 3 Smolts
  # 4 O1
  # 5 O2
  
  A <- matrix(0, 5, 5)
  colnames(A) <- rownames(A) <-
    c("Spawners","Parr","Smolts","O1","O2")
  
  # Smolts -> O1 via Grays Harbor
  A["O1","Smolts"] <- gh_surv
  
  # Ocean year 1
  A["Spawners","O1"] <- ocean_surv[1] * b2 * prespawn_surv
  A["O2","O1"]       <- ocean_surv[1] * (1 - b2)
  
  # Ocean year 2 (age-3 adults return)
  A["Spawners","O2"] <- ocean_surv[2] * prespawn_surv
  
  A
}

step_coho <- function(
    state,
    p_spawn_parr, c_spawn_parr,
    p_parr_smolt, c_parr_smolt,
    A_marine
) {
  
  # Freshwater DD stages
  parr_next <- beverton_holt(state["Spawners"],
                             p_spawn_parr,
                             c_spawn_parr)
  
  smolt_next <- beverton_holt(state["Parr"],
                              p_parr_smolt,
                              c_parr_smolt)
  
  # Marine linear transitions
  next_state <- as.numeric(A_marine %*% state)
  names(next_state) <- rownames(A_marine)
  
  # overwrite freshwater DD outputs
  next_state["Parr"]   <- parr_next
  next_state["Smolts"] <- smolt_next
  
  next_state
}

simulate_coho <- function(
    nyears,
    init = c(Spawners=500, Parr=0, Smolts=0, O1=0, O2=0),
    
    # Freshwater DD parameters (from EDT)
    p_spawn_parr,
    c_spawn_parr,
    p_parr_smolt,
    c_parr_smolt,
    
    # Marine parameters
    SAR,
    prespawn_surv
) {
  
  A <- make_A_marine_coho(
    SAR = SAR,
    prespawn_surv = prespawn_surv
  )
  
  out <- matrix(NA, nyears+1, 5)
  colnames(out) <- names(init)
  out[1,] <- init
  
  state <- init
  for(t in 1:nyears){
    state <- step_coho(
      state,
      p_spawn_parr,
      c_spawn_parr,
      p_parr_smolt,
      c_parr_smolt,
      A
    )
    out[t+1,] <- state
  }
  
  as.data.frame(out) |>
    transform(Year = 0:nyears, .before = 1)
}


#---- test with edt values 

for(sp in subpops){
  for(sc in scenarios){
    
    dat_sub <- dat %>%
      filter(Subpopulation == sp,
             Scenario == sc)
    
    if(nrow(dat_sub) == 0) next

Eggs.p.spawner <- 1500
p_spawn_parr  <- dat_sub$Productivity
c_spawn_parr  <- dat_sub$Capacity

p_parr_smolt  <- dat_sub$`Juvenile Productivity`
c_parr_smolt  <- dat_sub$`Juvenile Capacity`

SAR           <- dat_sub$`Juvenile-to-Adult Survival`

prespawn_surv <- 0.87505   # constant from your code

sim_coho <- simulate_coho(
  nyears = 100,
  
  p_spawn_parr  = p_spawn_parr,
  c_spawn_parr  = c_spawn_parr,
  
  p_parr_smolt  = p_parr_smolt,
  c_parr_smolt  = c_parr_smolt,
  
  SAR = SAR,
  prespawn_surv = prespawn_surv
)

#Compute Equilibrium
eq_mu <- sim_coho %>%
  filter(Year > 80) %>%
  summarise(eq_mu = mean(Spawners)) %>%
  pull(eq_mu)

#Store Results
results_eq_coho <- rbind(
  results_eq_coho,
  data.frame(
    Subpopulation = sp,
    Scenario = sc,
    Equilibrium_Mu_Spawners = eq_mu
  )
)
  }
}


eq_input <- dat %>%
  select(Subpopulation, Scenario, `Equilibrium Abundance`) %>%
  distinct()

results_eq_coho <- results_eq_coho %>%
  left_join(eq_input, by = c("Subpopulation","Scenario")) %>%
  mutate(
    Percent_Difference =
      100 * (Equilibrium_Mu_Spawners - `Equilibrium Abundance`) /
      `Equilibrium Abundance`
  )


