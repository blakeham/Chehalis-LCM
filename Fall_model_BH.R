# ============================================================
# FALL CHINOOK EDT → LIFE CYCLE MODEL (LCM) BATCH RUNNER
# ============================================================
#
# PURPOSE:
#   For each Fall Chinook row in the EDT dataframe, run a
#   Beverton-Holt + matrix marine life cycle model for 100 years
#   and return the equilibrium adult spawner abundance.
#
# MODEL STRUCTURE:
#   Freshwater: Density-Dependent Beverton-Holt recruitment
#               Spawners → Smolts
#
#   Marine:     Density-Independent linear transition matrix
#               Smolts → Ocean Years → Adult Returns
#
# KEY ASSUMPTIONS:
#   - Recruitment nonlinearity occurs ONLY at spawning.
#   - Marine survival is linear and age-structured.
#   - Fish that do not mature at age-5 exit the model.
#   - 100 years is assumed sufficient for equilibrium.
#
# ============================================================

edt <- read.csv("C:/Users/Blake.Hamilton/OneDrive - Kleinschmidt Associates/Projects/Chehalis/Pop ass model/EDT Outputs/Master EDT outputs.csv")


# =========================
# RECRUITMENT FUNCTION
# =========================
# Beverton-Holt density-dependent recruitment.
# N = number of spawners
# p = juvenile productivity
# c = juvenile capacity
#
# Produces next-year smolts.

beverton_holt <- function(N, p, c) {
  (p * N) / (1 + (p / c) * N)
}



# =========================
# MARINE TRANSITION MATRIX
# =========================
# Builds the linear density-independent transition matrix.
#
# State vector order:
#   1 Spawners
#   2 Smolts
#   3 O1 (Ocean year 1)
#   4 O2
#   5 O3
#   6 O4
#   7 O5
#
# SAR = smolt-to-adult return from EDT
# prespawn_surv = adult survival before spawning
#
# ocean_surv = annual marine survival
# b = maturation probabilities at ages 2-5

make_A_marine_fall_chinook <- function(
    SAR,
    prespawn_surv,
    ocean_surv = c(0.6, 0.7, 0.8, 0.9, 0.9),
    b = c(b2 = 0.007, b3 = 0.07, b4 = 0.35, b5 = 0.85)
) {
  
  # Back-calculate estuary/Grays Harbor survival
  gh_surv <- SAR / prod(ocean_surv)
  
  # Initialize empty matrix
  A <- matrix(0, nrow = 7, ncol = 7)
  colnames(A) <- rownames(A) <- c("Spawners","Smolts","O1","O2","O3","O4","O5")
  
  # Smolts → Ocean Year 1
  A["O1","Smolts"] <- gh_surv
  
  # Age-2 return split
  A["Spawners","O1"] <- ocean_surv[1] * b["b2"] * prespawn_surv
  A["O2","O1"]       <- ocean_surv[1] * (1 - b["b2"])
  
  # Age-3 return split
  A["Spawners","O2"] <- ocean_surv[2] * b["b3"] * prespawn_surv
  A["O3","O2"]       <- ocean_surv[2] * (1 - b["b3"])
  
  # Age-4 return split
  A["Spawners","O3"] <- ocean_surv[3] * b["b4"] * prespawn_surv
  A["O4","O3"]       <- ocean_surv[3] * (1 - b["b4"])
  
  # Age-5 return split
  A["Spawners","O4"] <- ocean_surv[4] * b["b5"] * prespawn_surv
  A["O5","O4"]       <- ocean_surv[4] * (1 - b["b5"])
  
  # Age-6+ fish exit the system
  A
}



# =========================
# ONE TIME STEP
# =========================
# Performs one yearly update.

step_fall_chinook <- function(state, p_spawn_smolt, c_spawn_smolt, A_marine) {
  
  # Density-dependent recruitment
  smolts_next <- beverton_holt(state["Spawners"], p_spawn_smolt, c_spawn_smolt)
  
  # Insert smolts into marine projection
  state_for_A <- state
  state_for_A["Smolts"] <- smolts_next
  
  # Linear marine transitions
  next_state <- as.numeric(A_marine %*% state_for_A)
  names(next_state) <- rownames(A_marine)
  
  # Keep smolt count for reporting
  next_state["Smolts"] <- smolts_next
  
  next_state
}



# =========================
# FULL SIMULATION
# =========================
# Runs model forward nyears and returns final spawners.

simulate_fall_chinook <- function(
    nyears,
    init,
    p_spawn_smolt,
    c_spawn_smolt,
    SAR,
    prespawn_surv
) {
  
  A <- make_A_marine_fall_chinook(SAR, prespawn_surv)
  
  state <- init
  for (t in 1:nyears) {
    state <- step_fall_chinook(state, p_spawn_smolt, c_spawn_smolt, A)
  }
  
  # Return equilibrium adult abundance
  state["Spawners"]
}



# =========================
# ROW-LEVEL RUNNER
# =========================
# Pulls EDT values and executes simulation.

run_fall_row <- function(row) {
  
  p_spawn_smolt <- as.numeric(row["Juvenile.Productivity"])
  c_spawn_smolt <- as.numeric(row["Juvenile.Capacity"])
  Juv_spawn_surv <- as.numeric(row["Juvenile.to.Spawner.Survival"])
  P  <- as.numeric(row["Productivity"])
  p1 <- as.numeric(row["Juvenile.to.Adult.Productivity"])
  
  # Data validity checks
  if (any(is.na(c(p_spawn_smolt, c_spawn_smolt, Juv_spawn_surv, P, p1)))) return(NA_real_)
  if (p1 == 0) return(NA_real_)
  
  # Pre-spawn survival
  prespawn_surv <- P / p1
  if (prespawn_surv <= 0) return(NA_real_)
  
  # Smolt-to-adult return ratio
  SAR <- Juv_spawn_surv / prespawn_surv
  if (SAR <= 0) return(NA_real_)
  
  simulate_fall_chinook(
    nyears = 100,
    init = c(Spawners = 200, Smolts = 0, O1 = 0, O2 = 0, O3 = 0, O4 = 0, O5 = 0),
    p_spawn_smolt = p_spawn_smolt,
    c_spawn_smolt = c_spawn_smolt,
    SAR = SAR,
    prespawn_surv = prespawn_surv
  )
}



# =========================
# EXECUTION BLOCK
# =========================

# 1) Create a NEW dataframe containing only Fall Chinook
edt_fall <- edt[edt$Species == "Fall Chinook", ]

# 2) Create output column in the NEW dataframe
edt_fall$Model_Equilibrium_Spawners <- NA_real_

# 3) Run model row-wise on edt_fall only
edt_fall$Model_Equilibrium_Spawners <-
  apply(edt_fall, 1, run_fall_row)

# 4) edt_fall now contains only Fall Chinook + model results
#    Original edt dataframe remains unchanged

# Optional: View quick summary
summary(edt_fall$Model_Equilibrium_Spawners)

# Optional: write to CSV
write.csv(edt_fall, "fall_chinook_model_results.csv", row.names = FALSE)

# ==============================
# Visualize Basin Trajectory 
# ==============================

simulate_fall_chinook_trajectory <- function(
    nyears,
    init,
    p_spawn_smolt,
    c_spawn_smolt,
    SAR,
    prespawn_surv
) {
  
  A <- make_A_marine_fall_chinook(SAR, prespawn_surv)
  
  # store yearly outputs
  out <- numeric(nyears + 1)
  out[1] <- init["Spawners"]
  
  state <- init
  for (t in 1:nyears) {
    state <- step_fall_chinook(state, p_spawn_smolt, c_spawn_smolt, A)
    out[t + 1] <- state["Spawners"]
  }
  
  data.frame(
    Year = 0:nyears,
    Spawners = out
  )
}

plot_rows <- edt[
  edt$Species == "Fall Chinook" &
    edt$Subpopulation == "ALL SUBPOPULATIONS",
]

traj_list <- list()

for (i in 1:nrow(plot_rows)) {
  
  row <- plot_rows[i, ]
  
  p_spawn_smolt <- row$Juvenile.Productivity
  c_spawn_smolt <- row$Juvenile.Capacity
  Juv_spawn_surv <- row$Juvenile.to.Spawner.Survival
  P  <- row$Productivity
  p1 <- row$Juvenile.to.Adult.Productivity
  
  prespawn_surv <- P / p1
  SAR <- Juv_spawn_surv / prespawn_surv
  
  traj <- simulate_fall_chinook_trajectory(
    nyears = 100,
    init = c(Spawners = 200, Smolts = 0, O1 = 0, O2 = 0, O3 = 0, O4 = 0, O5 = 0),
    p_spawn_smolt = p_spawn_smolt,
    c_spawn_smolt = c_spawn_smolt,
    SAR = SAR,
    prespawn_surv = prespawn_surv
  )
  
  traj$Scenario <- row$Scenario
  traj_list[[i]] <- traj
}

traj_df <- do.call(rbind, traj_list)


#Plot
library(ggplot2)

ggplot(traj_df, aes(Year, Spawners, color = Scenario)) +
  geom_line(size = 1.2) +
  theme_bw() +
  labs(
    title = "Fall Chinook – ALL SUBPOPULATIONS Trajectories",
    y = "Adult Spawners"
  )
