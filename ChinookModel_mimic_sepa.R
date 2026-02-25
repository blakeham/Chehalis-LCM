#===================================
#
#   Chinook model for SEPA
#
#
#
#
#===================================




# Build a linear transition matrix for the DI (marine) part.
# State vector order:
#   1 Spawners
#   2 Smolts
#   3 O1  (ocean year 1 fish; "Age 1" box in Fig 2)
#   4 O2  (ocean year 2 fish; potential age-2 returns)
#   5 O3  (ocean year 3 fish; potential age-3 returns)
#   6 O4  (ocean year 4 fish; potential age-4 returns)
#   7 O5  (ocean year 5 fish; potential age-5 returns)
#
# NOTE: Recruitment (Spawners -> Smolts) is nonlinear DD, so we keep that outside A.

#do no vary marine survival or maturity rates
make_A_marine_fall_chinook <- function(
    SAR,                   # total smolt-to-adult return from EDT (DI)
    prespawn_surv,          # pre-spawn survival p (DI)
    ocean_surv = c(0.6, 0.7, 0.8, 0.9, 0.9),  # Table 2 Fall Chinook
    b = c(b2 = 0.007, b3 = 0.07, b4 = 0.35, b5 = 0.85) # Table 2 Fall Chinook
) {
  stopifnot(length(ocean_surv) == 5)
  stopifnot(all(names(b) %in% c("b2","b3","b4","b5")))
  
  # Per report QA/QC description: Grays Harbor survival partitioned from SAR by
  # dividing SAR by product of annual ocean survivals (page 24).
  #maturity weighted survival
  b2 <- b["b2"]; b3 <- b["b3"]; b4 <- b["b4"]; b5 <- b["b5"]
  o1 <- ocean_surv[1]; o2 <- ocean_surv[2]; o3 <- ocean_surv[3]; o4 <- ocean_surv[4]
  
  W <- (o1 * b2) +
    (o1 * (1 - b2) * o2 * b3) +
    (o1 * (1 - b2) * o2 * (1 - b3) * o3 * b4) +
    (o1 * (1 - b2) * o2 * (1 - b3) * o3 * (1 - b4) * o4 * b5)
  
  gh_surv <- SAR / W
  
  #gh_surv <- SAR / prod(ocean_surv)
  
  # 7x7 matrix for DI transitions (we will overwrite nonlinear parts each step)
  A <- matrix(0, nrow = 7, ncol = 7)
  colnames(A) <- rownames(A) <- c("Spawners","Smolts","O1","O2","O3","O4","O5")
  
  # Smolts -> O1 via Grays Harbor survival (DI)
  A["O1", "Smolts"] <- gh_surv
  
  # Ocean year 1 survival: O1 -> O2, then maturity at age-2 splits to spawners vs remain in O2
  # Implement split after survival by putting terms into target rows:
  A["Spawners","O1"] <- ocean_surv[1] * b["b2"] * prespawn_surv
  A["O2",      "O1"] <- ocean_surv[1] * (1 - b["b2"])
  
  # Ocean year 2 survival: O2 -> O3, split by b3
  A["Spawners","O2"] <- ocean_surv[2] * b["b3"] * prespawn_surv
  A["O3",      "O2"] <- ocean_surv[2] * (1 - b["b3"])
  
  # Ocean year 3 survival: O3 -> O4, split by b4
  A["Spawners","O3"] <- ocean_surv[3] * b["b4"] * prespawn_surv
  A["O4",      "O3"] <- ocean_surv[3] * (1 - b["b4"])
  
  # Ocean year 4 survival: O4 -> O5, split by b5
  A["Spawners","O4"] <- ocean_surv[4] * b["b5"] * prespawn_surv
  A["O5",      "O4"] <- ocean_surv[4] * (1 - b["b5"])
  
  # Ocean year 5 survival exists in Table 2; without a b6 in the report,
  # we treat age-5+ non-maturing fish as exiting (no further state).
  # If you want them to die after year 5, do nothing; if you want them to return at age 6,
  # you'd need an additional maturation parameter (not provided in Table 2).
  # (So we do not advance O5.)
  A
}

#need to loop through prod and cap values as determined by draw
# 2 yr flood -> 1-p10-p100

step_fall_chinook <- function(
    state,
    p_spawn_smolt, c_spawn_smolt,
    A_marine
) {
  stopifnot(length(state) == 7)
  
  # 1) density-dependent recruitment (new smolts)
  smolts_next <- beverton_holt(state["Spawners"], p = p_spawn_smolt, c = c_spawn_smolt)
  
  # 2) IMPORTANT: let those new smolts enter the marine projection this time step
  state_for_A <- state
  state_for_A["Smolts"] <- smolts_next
  
  # 3) linear marine + returns + aging
  next_state <- as.numeric(A_marine %*% state_for_A)
  names(next_state) <- rownames(A_marine)
  
  # 4) keep Smolts as the newly produced smolts for reporting
  next_state["Smolts"] <- smolts_next
  
  next_state
}


#===== base sim function ===================
simulate_fall_chinook <- function(
    nyears,
    init = c(Spawners = 500, Smolts = 0, O1 = 0, O2 = 0, O3 = 0, O4 = 0, O5 = 0),
    # Freshwater DD parameters (from EDT; placeholders here)
    p_spawn_smolt,
    c_spawn_smolt,
    # Marine/pre-spawn DI parameters
    SAR,
    prespawn_surv,
    ocean_surv = c(0.6, 0.7, 0.8, 0.9, 0.9),
    b = c(b2 = 0.007, b3 = 0.07, b4 = 0.35, b5 = 0.85)
) {
  A <- make_A_marine_fall_chinook(
    SAR = SAR,
    prespawn_surv = prespawn_surv,
    ocean_surv = ocean_surv,
    b = b
  )
  
  out <- matrix(NA_real_, nrow = nyears + 1, ncol = 7)
  colnames(out) <- names(init)
  out[1, ] <- init
  state <- init
  
  for (t in 1:nyears) {
    state <- step_fall_chinook(
      state = state,
      p_spawn_smolt = p_spawn_smolt[t],
      c_spawn_smolt = c_spawn_smolt[t],
      A_marine = A
    )
    out[t + 1, ] <- state
  }
  
  as.data.frame(out) |>
    transform(Year = 0:nyears, .before = 1)
}
#=============================================================


# loop sim function
simulate_fall_chinook_tv_const_ocean <- function(
    nyears,
    init = c(Spawners=500, Smolts=0, O1=0, O2=0, O3=0, O4=0, O5=0),
    params,  # data.frame with columns: p_spawn_smolt, c_spawn_smolt, prespawn_surv, Juv_spawn_surv
    
    ocean_surv = c(0.6, 0.7, 0.8, 0.9, 0.9),  # constant
    b = c(b2 = 0.007, b3 = 0.07, b4 = 0.35, b5 = 0.85) # constant
) {
  stopifnot(nrow(params) >= nyears)
  
  out <- matrix(NA_real_, nrow = nyears + 1, ncol = length(init))
  colnames(out) <- names(init)
  out[1,] <- init
  state <- init
  
  for (t in 1:nyears) {
    
    # 1) year-specific freshwater DD recruitment
    smolts_next <- beverton_holt(
      state["Spawners"],
      p = params$p_spawn_smolt[t],
      c = params$c_spawn_smolt[t]
    )
    
    # 2) year-specific SAR (pre-ocean) / prespawn; ocean survivals fixed
    prespawn_t <- params$prespawn_surv[t]
    SAR_t <- params$Juv_spawn_surv[t] / prespawn_t
    
    A_t <- make_A_marine_fall_chinook(
      SAR = SAR_t,
      prespawn_surv = prespawn_t,
      ocean_surv = ocean_surv,
      b = b
    )
    
    # 3) advance marine states
    state_for_A <- state
    state_for_A["Smolts"] <- smolts_next
    next_state <- as.numeric(A_t %*% state_for_A)
    names(next_state) <- rownames(A_t)
    
    # 4) keep smolts for reporting
    next_state["Smolts"] <- smolts_next
    
    state <- next_state
    out[t+1,] <- state
  }
  
  as.data.frame(out) |>
    transform(Year = 0:nyears, .before = 1)
}

#=================== end loop sim function ======================


# ---- Example run (you MUST replace these EDT-derived placeholders) ----
# p_spawn_smolt and c_spawn_smolt should come from EDT for the scenario/flood type you want.
# SAR is EDT "Total returns to Spawners" DI ratio (Table 2 descriptor), and prespawn_surv is EDT p.

#pull values from dataframe
#p_spawn_smolt -> # Juvenile Productivity
#c_spawn_smolt -> # Juvenile Capacity
#prespawn_surv -> # Prespawner Survival
#Juv_spawn_surv -> # Juvenile-to-Spawner Survival this column has a percent on it

#to do draws could create vector of values that align with t draws
# Using your provided row values
##point to sepa for test run
dat <- read_xlsx('C:/Users/Jason.Romine/OneDrive - Kleinschmidt Associates/Chehalis_LCM_devo/LCM Input and Results files/Chinook/data/All_Species_Median_Chinook.xlsx', 
                 sheet = 'Values')

dat <- read.csv("C:/Users/Blake.Hamilton/OneDrive - Kleinschmidt Associates/Projects/Chehalis/Pop ass model/SEPA LCM/Chinook/data/Master SEPA EDT Outputs.csv",
                check.names = FALSE)
#we only care about above crim and rbf to crim fall chinook and near for this
# dat <- dat %>% filter(Subpopulation %in% c('Above Crim SB', 'Chehalis RBF to Crim SB'),
#                       Species == 'Fall Chinook',
#                       row_number() %in% grep("Near", Scenario, ignore.case = T),
#                       row_number() %in% grep('NA', Scenario, ignore.case = T))


#above crim near, No Action
dat <- dat %>%
  filter(
    Subpopulation %in% c("Above Crim SB", "Chehalis RBF to Crim SB"),
    Species == "Fall Chinook"
  )

#Get unique Scenarios
scenarios <- unique(dat$Scenario)
subpops <- c("Above Crim SB", "Chehalis RBF to Crim SB")

results_eq <- data.frame()

#now we have 2, 10 and 100 yr flood
#for near there really isn't much difference in values
#lets just sim the near 2yr over time and ess what happens

#Make a flood scenario
flood <- c('2yr', '10yr', '100yr')

# flood_prob_table <- tibble::tribble(
#   ~period,                  ~typical, ~major, ~catastrophic,
#   "recent",                    0.85,    0.14,     0.01,
#   "ensemble_median_mid",       0.71,    0.27,     0.016,
#   "ensemble_median_late",      0.63,    0.34,     0.027,
#   "ensemble_max_mid",          0.58,    0.37,     0.05,
#   "ensemble_max_late",         0.45,    0.47,     0.076
# )

#for sepa comps use '2yr', '10yr', '100yr'
flood_prob_table <- tibble::tribble(
  ~period,                  ~'2yr', ~'10yr', ~'100yr',
  "recent",                    0.85,    0.14,     0.01,
  "ensemble_median_mid",       0.71,    0.27,     0.016,
  "ensemble_median_late",      0.63,    0.34,     0.027,
  "ensemble_max_mid",          0.58,    0.37,     0.05,
  "ensemble_max_late",         0.45,    0.47,     0.076
)



# quick check: rows sum to 1 (allowing tiny rounding error)
#stopifnot(all(abs(rowSums(flood_prob_table[,c("typical","major","catastrophic")]) - 1) < 1e-6))

#note the recent or other here must change for future
draw_flood_types <- function(nyears, period = "recent", prob_table = flood_prob_table, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  row <- prob_table[prob_table$period == period, ]
  if (nrow(row) != 1) stop("period not found in prob_table: ", period)
  
  probs <- as.numeric(row[, c("2yr","10yr","100yr")])
  types <- c("2yr","10yr","100yr")
  
  sample(types, size = nyears, replace = TRUE, prob = probs)
}



#------------------------------- loop me ---------------------------
for(sp in subpops){
  for(sc in scenarios){
    
    dat_sub <- dat %>%
      filter(Subpopulation == sp,
             Scenario == sc)
    
    if(nrow(dat_sub) == 0) next

for(j in 1:500)
{

# Draw flood types:
#set the time frame for the run here
nyears = 75

flood_seq <- draw_flood_types(nyears = nyears, period = "recent")
table(flood_seq) / length(flood_seq)
flood_seq

p_spawn_smolt = rep(NA, nyears)
c_spawn_smolt = rep(NA, nyears)
prespawn_surv = rep(NA, nyears) 
Juv_spawn_surv = rep(NA, nyears)

for(i in 1:nyears)
{
  p_spawn_smolt[i]  <- dat_sub$`Juvenile Productivity`
  c_spawn_smolt[i]  <- dat_sub$`Juvenile Capacity`
  prespawn_surv[i]  <- dat_sub$`Prespawner Survival`
  Juv_spawn_surv[i] <- dat_sub$`Juvenile-to-Spawner Survival`
}

#set params here for function to step through; might be more efficient to have the function draw them

params <- data.frame(Year = 1:nyears, 
                     p_spawn_smolt,
                     c_spawn_smolt,
                     prespawn_surv, 
                     Juv_spawn_surv,
                     o1=0.6, o2=0.7, o3=0.8, o4=0.9, o5=0.9,
                     b2=0.007, b3=0.07, b4=0.35, b5=0.85)


# Make SAR consistent with Juvenile-to-Spawner survival = 0.027
# so that SAR * prespawn_surv ≈ 0.027
SAR <- Juv_spawn_surv / prespawn_surv  # ≈ 0.03242


#need to capture each sim into aframe
sim <- simulate_fall_chinook_tv_const_ocean(
  nyears = nyears,
  init = c(Spawners = 310, Smolts = 15000, O1 = 1300, O2 = 775, O3 = 505, O4 = 262, O5 = 35),
  params = params
  
)

#capture the simulations
sim$simrun = j
if(j == 1){all_sims = sim}else{all_sims =  rbind(all_sims, sim)}


#tail(sim)
#plot the first run then do line
if(j == 1){plot(sim$Year, sim$Spawners, type="l", xlab="Year", ylab="Spawners", ylim = c(300,325))}else{
lines(sim$Year, sim$Spawners, col = 'azure4')}

#add avg line

}

#add avg line

mu_spawners <- all_sims %>%
  group_by(Year) %>%
  summarise(mu = mean(Spawners))

with(mu_spawners, lines(Year, mu, lwd = 2, col = 1))

#Define equilibrium as mean of last 20 years of sim
eq_mu <- mu_spawners %>%
  filter(Year > (nyears - 20)) %>%
  summarise(eq_mu = mean(mu)) %>%
  pull(eq_mu)

#Put all equilibrium values for each scenario/subpop into dataframe
results_eq <- rbind(
  results_eq,
  data.frame(
    Subpopulation = sp,
    Scenario = sc,
    Equilibrium_Mu_Spawners = eq_mu
  )
)
  }
}

#Pull EDT input Equilibrium Abundance value to compare to model output
eq_input <- dat %>%
  select(Subpopulation, Scenario, `Equilibrium Abundance`) %>%
  distinct()

#Add in the Equilibrium abundance from dat and compute % diff
results_eq <- results_eq %>%
  left_join(eq_input, by = c("Subpopulation", "Scenario")) %>%
  mutate(
    Percent_Difference = 
      100 * (Equilibrium_Mu_Spawners - `Equilibrium Abundance`) /
      `Equilibrium Abundance`
  )
