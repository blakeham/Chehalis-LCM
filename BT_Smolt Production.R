# ==========================================================
# Fall Chinook LCM: Adults → Smolts to Grays Harbor
# FRE Operations vs No Action
#
#
# This model was developed to be as simple as possible and provide
# liberal estimates of mortality associated with FRE operations to 
# quantify the total difference in smolts reaching Grays Harbor
#
# ==========================================================

library(ggplot2)

# ------------------------------
# 1. Initial Conditions
# ------------------------------

#Want to only apply the flood-based, FRE ops mortalities to the proportion of the
#population above Newaukum (minus S Fork Chehalis), which we are calling Upper Chehalis
#Entire basin pop based on ASRP TM2 = 29,358 fish
#Upper Chehalis Basin pop based on ASRP TM2: 
          #SF to Elk(187) + Elk(19) + Elk to Crim(85) + Above Crim(113) = 404 fish under current conditions
#All other areas (S Fork + all downstream) = 28,954
#So, the area of focus is only 1.40% of all adults in system

pop_init <- 29358  
n_years  <- 70
n_mc     <- 500
set.seed(123)

# ------------------------------
# 2. Life-History Parameters
# ------------------------------
adults_per_redd     <- 2.5 #Washington State Salmon Recovery Technical Memoranda
eggs_per_redd     <- 5400 #NOAA LCM

#
egg_to_fry        <- 0.30   # Egg-to-fry survival for Chinook salmon typically ranges from 25–70% depending on incubation conditions,
                            #with life-cycle models commonly assuming values around 0.30–0.40. This model uses 0.40 as a
                            #representative mid-range estimate consistent with NOAA life-cycle modeling assumptions.
fry_to_smolt_fw   <- 0.15   # Can pull older info from: https://www.waterboards.ca.gov/waterrights/water_issues/programs/bay_delta/deltaflow/docs/exhibits/nmfs/spprt_docs/nmfs_exh4_healey_1991.pdf
outmigrant_surv      <- 0.70   # 

egg_capacity      <- 8e6

BH <- function(S, p, c){
  (S * p) / (1 + (p / c) * S)
}

# ------------------------------
# 3. Flood Probability Regime
# ------------------------------

#Use ensemble maximum values for most conservative estimate - higher mortality
flood_probs <- data.frame(
  period = c("Mid","Late"),
  typical = c(0.58,0.45),
  major = c(0.37,0.47),
  catastrophic = c(0.05,0.076)
)

get_climate_period <- function(year){
  if(year <= 35) return("Mid")
  return("Late")
}

draw_flood <- function(year){
  period <- get_climate_period(year)
  probs <- flood_probs[flood_probs$period==period,
                       c("typical","major","catastrophic")]
  sample(c("typical","major","catastrophic"),
         1,
         prob=as.numeric(probs))
}

# ------------------------------
# 4. Scenario-Specific Redd Mortalities
# ------------------------------

#the 3.5% for major and 10.5% for catastrophic flood redd mortalities only apply
#to the portion of the pop above Newaukum (minus S Fork Chehalis). That is only 1.40% 
#of the total pop (404/29,358 * 100).

# 404 fish / 2.5 fish per redd = 161.6 redds in Upper Chehalis
# 10.5% of these wiped out with catastrophic flood: 161.6 * 0.105 = 17 redds lost
# 3.5% of these wiped out with catastrophic flood: 161.6 * 0.035 = 6 redds lost

#Assume total Chehalis pop redds is 29,358 / 2.5 fish per redd = 11,743 redds
#Upper Chehalis redd loss as proportion of entire basin: 
            #Catastrophic flood: 17 redds lost / 11,743 total redds = 0.00145
            #Major flood: 6 redds lost / 11,743 = 0.0005
#Apply these proportions for the redd mort with FRE ops

redd_mort_dam <- c(
  typical = 0.00,
  major = 0.0005,
  catastrophic = 0.00145
)

redd_mort_no_action <- c(
  typical = 0.00,
  major = 0.000,
  catastrophic = 0.000
)

# ------------------------------
# 5. Simulation Function
# ------------------------------
run_lcm <- function(redd_mort_vector){
  
  smolts_out <- matrix(NA,n_mc,n_years)
  
  for(mc in 1:n_mc){
    
    adults <- pop_init
    
    for(yr in 1:n_years){
      
      redds   <- adults / adults_per_redd
      
      flood_type <- draw_flood(yr)
      redds_post <- redds * (1 - redd_mort_vector[flood_type])
      
      eggs <- BH(redds_post, eggs_per_redd, egg_capacity)
      fry  <- eggs * egg_to_fry
      smolts_fw <- fry * fry_to_smolt_fw
      smolts_GH <- smolts_fw * outmigrant_surv
      
      smolts_out[mc,yr] <- smolts_GH
      adults <- pop_init
    }
  }
  
  smolts_out
}

# ------------------------------
# 6. Run Both Scenarios
# ------------------------------
results_dam       <- run_lcm(redd_mort_dam)
results_no_action <- run_lcm(redd_mort_no_action)

# ------------------------------
# 7. Percent Decline Calculation
# ------------------------------
mean_dam       <- apply(results_dam,2,mean)
mean_no_action <- apply(results_no_action,2,mean)

percent_decline_by_year <- 
  (mean_no_action - mean_dam) / mean_no_action * 100

mean_percent_decline <- mean(percent_decline_by_year)

cat("Mean Percent Decline in Smolts (FRE vs No Action):",
    round(mean_percent_decline,2), "%\n")

# ------------------------------
# 7b. Summary Comparison Table
# ------------------------------

# Mean smolts across all years (and MC runs)
overall_mean_dam       <- mean(results_dam)
overall_mean_no_action <- mean(results_no_action)

# Absolute and percent difference
difference_smolt <- overall_mean_no_action - overall_mean_dam
percent_diff     <- difference_smolt / overall_mean_no_action * 100

summary_table <- data.frame(
  Scenario = c("No Action", "FRE Operations"),
  Mean_Smolts = c(overall_mean_no_action,
                  overall_mean_dam)
)

comparison_row <- data.frame(
  Scenario = "Difference (No Action - FRE)",
  Mean_Smolts = difference_smolt
)

percent_row <- data.frame(
  Scenario = "Percent Difference (%)",
  Mean_Smolts = percent_diff
)

summary_output <- rbind(summary_table,
                        comparison_row,
                        percent_row)

summary_output$Mean_Smolts <- round(summary_output$Mean_Smolts, 3)
print(summary_output)

# ------------------------------
# 8. Prepare Plot Data
# ------------------------------
long_dam <- as.data.frame(as.table(results_dam))
colnames(long_dam) <- c("Sim","Year","Smolts")
long_dam$Scenario <- "FRE Operations"

long_no_action <- as.data.frame(as.table(results_no_action))
colnames(long_no_action) <- c("Sim","Year","Smolts")
long_no_action$Scenario <- "No Action"

long_all <- rbind(long_dam, long_no_action)
long_all$Year <- as.numeric(gsub("V","", long_all$Year))

summary_df <- data.frame(
  Year = rep(1:n_years,2),
  Smolts = c(mean_dam, mean_no_action),
  Scenario = rep(c("FRE Operations","No Action"),
                 each=n_years)
)

# ------------------------------
# 9. Spaghetti Plot
# ------------------------------

ggplot() +
  
  # Mean lines (bold)
  geom_line(data=summary_df,
            aes(x=Year,y=Smolts,color=Scenario),
            linewidth=1.4) +
  
  scale_color_manual(values=c(
    "FRE Operations"="red",
    "No Action"="blue"
  )) +
  
  labs(title="Fall Chinook Smolts Reaching Grays Harbor",
       subtitle=paste("Decline with FRE:",
                      round(comparison_row$Mean_Smolts,2),"smolts;",
                      round(summary_output$Mean_Smolts[4],4), "%"),
       x="Year",
       y="Smolts to Grays Harbor",
       color="Scenario") +
  
  theme_minimal(base_size=14)
