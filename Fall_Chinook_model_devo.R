#==================================
# 
# 
#  LCM based on SEPA LCM utilizes EDT outputs:
#  productivity
#  capacity
#  prespawner survival
#  juvenile to spawner survival
# 
# 
#
#=================================



#===  Functions  

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

g <- function(x) head(x)



# Fall Chinook matrix-type life cycle model (Figure 2 + Table 2)
# - Nonlinear DD Beverton-Holt for Spawners -> Smolts
# - DI (linear) survival through Grays Harbor + ocean years
# - Age-specific maturation b2-b5 producing spawning adults
# Source: Beechie & Hendrix 2025 SEPA EDT-LCM report (Table 2, Figure 2)
#  adults that do not mature at b5 are killed

beverton_holt <- function(N, p, c) {
  # N_{t+1} = (p*N) / (1 + (p/c)*N)
  (p * N) / (1 + (p / c) * N)
}

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
  gh_surv <- SAR / prod(ocean_surv)
  
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
      p_spawn_smolt = p_spawn_smolt,
      c_spawn_smolt = c_spawn_smolt,
      A_marine = A
    )
    out[t + 1, ] <- state
  }
  
  as.data.frame(out) |>
    transform(Year = 0:nyears, .before = 1)
}

# ---- Example run (you MUST replace these EDT-derived placeholders) ----
# p_spawn_smolt and c_spawn_smolt should come from EDT for the scenario/flood type you want.
# SAR is EDT "Total returns to Spawners" DI ratio (Table 2 descriptor), and prespawn_surv is EDT p.

#equilibirum/juvenile to adult abundance

#equilibrium abundance is column S
#juvenile to adult abundance is column AA
#SAR is equilibirum/juvenile to adult

# Using your provided row values
p_spawn_smolt <- 163.5      # Juvenile Productivity
c_spawn_smolt <- 5344       # Juvenile Capacity

prespawn_surv <- 0.832856728
Juv_spawn_surv <- 0.0270716299275952
# Make SAR consistent with Juvenile-to-Spawner survival = 0.027
# so that SAR * prespawn_surv ≈ 0.027
SAR <- Juv_spawn_surv / prespawn_surv  # ≈ 0.03242

sim <- simulate_fall_chinook(
  nyears = 80,
  init = c(Spawners = 200, Smolts = 0, O1 = 0, O2 = 0, O3 = 0, O4 = 0, O5 = 0),
  p_spawn_smolt = p_spawn_smolt,
  c_spawn_smolt = c_spawn_smolt,
  SAR = SAR,
  prespawn_surv = prespawn_surv
)

tail(sim)
plot(sim$Year, sim$Spawners, type="l", xlab="Year", ylab="Spawners")



#our edt

# Using your provided row values
p_spawn_smolt <- 223.1      # Juvenile Productivity
c_spawn_smolt <- 4541       # Juvenile Capacity

#must cals prespawn survial
#EDT out:
#spawner to spawer ->  productivty-> P = 4.5 column q column s in ours
P <-4.25161301642635
#juvi to adult productivity ->  p1 = 4.5 column Y
p1 <- 5.27369479553325

#p2 is pre-spawn survival 
p2 <- P/p1

prespawn_surv <- p2
#prespawn_surv <- 0.832856728
Juv_spawn_surv <- 0.0264235693250393
# Make SAR consistent with Juvenile-to-Spawner survival = 0.027
# so that SAR * prespawn_surv ≈ 0.027
SAR <- Juv_spawn_surv / prespawn_surv  # ≈ 0.03242

sim2 <- simulate_fall_chinook(
  nyears = 80,
  init = c(Spawners = 200, Smolts = 0, O1 = 0, O2 = 0, O3 = 0, O4 = 0, O5 = 0),
  p_spawn_smolt = p_spawn_smolt,
  c_spawn_smolt = c_spawn_smolt,
  SAR = SAR,
  prespawn_surv = prespawn_surv
)

tail(sim2)
lines(sim2$Year, sim2$Spawners, col = 2)

