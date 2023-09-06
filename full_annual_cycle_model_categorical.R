
library(tidyverse)
library(lubridate)

### prepare data ###
# read in data
spring_dat <- read.csv("spring_data.csv")
autumn_dat <- read.csv("autumn_data.csv")

# scale covariate data
spring_dat <- spring_dat %>% 
  mutate(scaled_mean_precip = as.numeric(scale(mean_precip)), 
         scaled_prop_days_below_freezing = as.numeric(scale(prop_days_below_freezing)), 
         scaled_prop_storm_days = as.numeric(scale(prop_storm_days)),
         scaled_prop_grass = as.numeric(scale(prop_grass)),
         scaled_prop_ag = as.numeric(scale(prop_ag)),
         scaled_prop_bog = as.numeric(scale(prop_bog))) %>% 
  # fill NAs for sub_seasons w/o these habitat types
  # (needed for model to run, but negated by fixing betas at 0)
  replace_na(list(scaled_prop_ag = 0, scaled_prop_grass = 0))

autumn_dat <- autumn_dat %>% 
  mutate(scaled_mean_precip = as.numeric(scale(mean_precip)), 
         scaled_prop_days_below_freezing = as.numeric(scale(prop_days_below_freezing)), 
         scaled_prop_storm_days = as.numeric(scale(prop_storm_days)), 
         scaled_prop_grass = as.numeric(scale(prop_grass)),
         scaled_prop_ag = as.numeric(scale(prop_ag)),
         scaled_prop_bog = as.numeric(scale(prop_bog))) %>% 
  # fill NAs for sub_seasons w/o these habitat types
  # (needed for model to run, but negated by fixing betas at 0)
  replace_na(list(scaled_prop_ag = 0, scaled_prop_grass = 0))

# subset data for indexing
# 1. a loop (spring sub-model, individual/year index)
spring_dat_annual <- spring_dat %>% filter(sub_season == 6)
spring_dat_annual <- spring_dat_annual %>% 
  rename(arrival_day = first_day) %>% 
  mutate(scaled_arrival_day = as.numeric(scale(arrival_day)))
# 2. b loop (spring model, individual/year/sub-season index)
spring_dat_all <- spring_dat
# 3. c loop (spring model, individual/year/(non migratory) sub-season index)
spring_dat_non_mig <- spring_dat %>% 
  filter(sub_season %in% c(1, 3, 4, 6))
# 4. d loop (autumn model, individual/year index)
autumn_dat_annual <- autumn_dat %>% 
  group_by(id, year) %>% 
  summarise(breeding_success = unique(breeding_success),
            breeding_outcome = unique(breeding_outcome),
            annual_survival = prod(survival), 
            max_sub_season = max(sub_season)) %>% 
  ungroup()
autumn_dat_annual <- merge(autumn_dat_annual,
                          (autumn_dat %>% 
                            filter(sub_season == 2) %>% 
                            rename(departure_day = first_day) %>% 
                            select(id, year, departure_day) %>% 
                            mutate(scaled_departure_day = as.numeric(scale(departure_day)))),
                         by = c("id", "year"), all.x = T)
# 5. e loop (autumn model, individual/year/sub-season index)
autumn_dat_all <- autumn_dat
# 6. f loop (autumn model, individual/year/(non migratory) sub-season index)
autumn_dat_non_mig <- autumn_dat %>% 
  filter(sub_season %in% c(1, 3, 4, 6))

### full annual cycle model in JAGS ###
# model code
cat(file = "full_annual_cycle_model_categorical.txt", "
model {
  
  #### priors
  
  ### spring sub-model
  ## sub seasons: 1 = late winter, 2 = migration flight from IE/GB to IS,
  ## 3 = 1/2 staging, 4 = 2/2 staging, 5 = migration flight from IS to GL,
  ## 6 = early breeding
  ## breeding outcomes: 1 = successful attempt, 2 = failed attempt,
  ## 3 = deferral
  
  ## effects of weather and habitat use on behavior
  # random intercepts for individual
  for(i in 1:nind){ # individual index
    alpha_ODBA_spring[i] ~ dnorm(0, tau_ODBA_spring)
    alpha_PTF_spring[i] ~ dnorm(0, tau_PTF_spring)
  } # i
  tau_ODBA_spring <- pow(sd_ODBA_spring, -2)
  sd_ODBA_spring ~ dunif(0, 1)
  tau_PTF_spring <- pow(sd_PTF_spring, -2)
  sd_PTF_spring ~ dunif(0, 1)
       
  # ODBA observation variation
  for(a in 1:nobs1){ # individual/year index
    tau_ODBA_obs_spring[ids1[a], years1[a]] <- pow(sd_ODBA_obs_spring[ids1[a], years1[a]], -2)
    sd_ODBA_obs_spring[ids1[a], years1[a]] ~ dunif(0, 10)
  } # a
  
  # betas
  for(s in 1:n_sub_seasons){ # sub-season index
    beta_ODBA_precip_spring[s] ~ dnorm(0, 1)
    beta_ODBA_temp_spring[s] ~ dnorm(0, 1)
  } # s
    
  for(s in c(1, 3, 4)){ # winter and staging sub-seasons
    beta_ODBA_storm_spring[s] <- 0
    beta_ODBA_grass_spring[s] ~ dnorm(0, 1)
    beta_ODBA_ag_spring[s] ~ dnorm(0, 1)
    beta_ODBA_bog_spring[s] ~ dnorm(0, 1)
    beta_PTF_precip_spring[s] ~ dnorm(0, 1)
    beta_PTF_temp_spring[s] ~ dnorm(0, 1)
    beta_PTF_grass_spring[s] ~ dnorm(0, 1)
    beta_PTF_ag_spring[s] ~ dnorm(0, 1)
    beta_PTF_bog_spring[s] ~ dnorm(0, 1)
  } # s
  
  for(s in c(2, 5)){ # migration sub-seasons
    beta_ODBA_storm_spring[s] ~ dnorm(0, 1)
    beta_ODBA_grass_spring[s] <- 0
    beta_ODBA_ag_spring[s] <- 0
    beta_ODBA_bog_spring[s] <- 0
    beta_PTF_precip_spring[s] <- 0
    beta_PTF_temp_spring[s] <- 0
    beta_PTF_grass_spring[s] <- 0
    beta_PTF_ag_spring[s] <- 0
    beta_PTF_bog_spring[s] <- 0
  } # s
  
  # breeding sub-season
  beta_ODBA_storm_spring[6] <- 0
  beta_ODBA_grass_spring[6] <- 0
  beta_ODBA_ag_spring[6] <- 0
  beta_ODBA_bog_spring[6] ~ dnorm(0, 1)
  beta_PTF_precip_spring[6] ~ dnorm(0, 1)
  beta_PTF_temp_spring[6] ~ dnorm(0, 1)
  beta_PTF_grass_spring[6] <- 0
  beta_PTF_ag_spring[6] <- 0
  beta_PTF_bog_spring[6] ~ dnorm(0, 1)
  
  ## effects of behavior on breeding outcome
  # random intercept for individual
  for(i in 1:nind){ # individual index
    for(k in 1:(n_B - 1)){ # breeding outcome index
      alpha_B_spring[i, k] ~ dnorm(0, tau_B_spring[k])
    } # k
    alpha_B_spring[i, n_B] <- 0
  } # i
  
  for(k in 1:(n_B - 1)){ # breeding outcome index
  # individual-state variation
    tau_B_spring[k] <- pow(sd_B_spring[k], -2)
    sd_B_spring[k] ~ dunif(0, 10)
  # mean betas
    mu_beta_ODBA_spring[k] ~ dnorm(0, 0.01)
    mu_beta_PTF_spring[k] ~ dnorm(0, 0.01)
  # sub-seasonal deviations for betas
    tau_beta_ODBA_spring[k] <- pow(sd_beta_ODBA_spring[k], -2)
    sd_beta_ODBA_spring[k] ~ dunif(0, 10)
    tau_beta_PTF_spring[k] <- pow(sd_beta_PTF_spring[k], -2)
    sd_beta_PTF_spring[k] ~ dunif(0, 10)
  # random variation for year effect
    tau_year_spring[k] <- pow(sd_year_spring[k], -2)
    sd_year_spring[k] ~ dunif(0, 10)
  # arrival date beta
    beta_arrival_breed[k] ~ dnorm(0, 0.01)
  } # k
  beta_arrival_breed[n_B] <- 0
  
  # sub-seasonal betas
  for(s in 1:n_sub_seasons){ # sub-season index
    for(k in 1:(n_B - 1)){ # breeding outcome index
      beta_ODBA_spring[s, k] ~ dnorm(mu_beta_ODBA_spring[k], tau_beta_ODBA_spring[k])
    } # k
    beta_ODBA_spring[s, n_B] <- 0
  } # s
  
  for(s in c(1, 3, 4, 6)){ # (non-migratory) sub-season index
    for(k in 1:(n_B - 1)){ # breeding outcome index
      beta_PTF_spring[s, k] ~ dnorm(mu_beta_PTF_spring[k], tau_beta_PTF_spring[k])
    } # k
    beta_PTF_spring[s, n_B] <- 0
  } # s
  
  # random year effect
  for(t in 1:nyears){ # year index
    for(k in 1:(n_B - 1)){ # breeding outcome index
      eps_year_spring[t, k] ~ dnorm(0, tau_year_spring[k])
    } # k
    eps_year_spring[t, n_B] <- 0
  } # t
  
  ### autumn sub-model
  ## sub seasons: 1 = late breeding, 2 = migration flight from GL to IS,
  ## 3 = 1/2 staging, 4 = 2/2 staging, 5 = migration flight from IS to IE/GB,
  ## 6 = early winter
  ## breeding outcomes: 1 = successful attempt, 2 = failed attempt,
  ## 3 = deferral
  
  ## effects of breeding outcome, weather, and habitat use on behavior
  # random intercepts for individual
  for(i in ind_autumn[]){ # individual index
    alpha_ODBA_autumn[i] ~ dnorm(0, tau_ODBA_autumn)
    alpha_PTF_autumn[i] ~ dnorm(0, tau_PTF_autumn)
  } # i
  tau_ODBA_autumn <- pow(sd_ODBA_autumn, -2)
  sd_ODBA_autumn ~ dunif(0, 1)
  tau_PTF_autumn <- pow(sd_PTF_autumn, -2)
  sd_PTF_autumn ~ dunif(0, 1)
  
  # observation variation
  for(d in 1:nobs4){ # individual/year index
    tau_ODBA_obs_autumn[ids4[d], years4[d]] <- pow(sd_ODBA_obs_autumn[ids4[d], years4[d]], -2)
    sd_ODBA_obs_autumn[ids4[d], years4[d]] ~ dunif(0, 10)
  } # d
  
  # betas
  for(s in 1:n_sub_seasons){ # sub-season index
    for(k in 2:n_B){ # breeding outcome index
      beta_ODBA_B_autumn[s, k] ~ dnorm(0, 1)
    } # k
    beta_ODBA_B_autumn[s, 1] <- 0
    beta_ODBA_precip_autumn[s] ~ dnorm(0, 1)
    beta_ODBA_temp_autumn[s] ~ dnorm(0, 1)
  } # s
  
  for(s in c(3, 4, 6)){ # staging and winter sub-seasons
    beta_ODBA_storm_autumn[s] <- 0
    beta_ODBA_grass_autumn[s] ~ dnorm(0, 1)
    beta_ODBA_ag_autumn[s] ~ dnorm(0, 1)
    beta_ODBA_bog_autumn[s] ~ dnorm(0, 1)
    for(k in 2:n_B){ # breeding outcome index
      beta_PTF_B_autumn[s, k] ~ dnorm(0, 1)
    } # k
    beta_PTF_B_autumn[s, 1] <- 0
    beta_PTF_precip_autumn[s] ~ dnorm(0, 1)
    beta_PTF_temp_autumn[s] ~ dnorm(0, 1)
    beta_PTF_grass_autumn[s] ~ dnorm(0, 1)
    beta_PTF_ag_autumn[s] ~ dnorm(0, 1)
    beta_PTF_bog_autumn[s] ~ dnorm(0, 1)
  } # s
  
  for(s in c(2, 5)){ # migratory sub-seasons
    beta_ODBA_storm_autumn[s] ~ dnorm(0, 1)
    beta_ODBA_grass_autumn[s] <- 0
    beta_ODBA_ag_autumn[s] <- 0
    beta_ODBA_bog_autumn[s] <- 0
    for(k in 1:n_B){ # breeding outcome index
      beta_PTF_B_autumn[s, k] <- 0
    } # k
    beta_PTF_precip_autumn[s] <- 0
    beta_PTF_temp_autumn[s] <- 0
    beta_PTF_grass_autumn[s] <- 0
    beta_PTF_ag_autumn[s] <- 0
    beta_PTF_bog_autumn[s] <- 0
  } # s
  
  # breeding sub-season
  beta_ODBA_storm_autumn[1] <- 0
  beta_ODBA_grass_autumn[1] <- 0
  beta_ODBA_ag_autumn[1] <- 0
  beta_ODBA_bog_autumn[1] ~ dnorm(0, 1)
  for(k in 2:n_B){ # breeding outcome index
    beta_PTF_B_autumn[1, k] ~ dnorm(0, 1)
  } # k
  beta_PTF_B_autumn[1, 1] <- 0
  beta_PTF_precip_autumn[1] ~ dnorm(0, 1)
  beta_PTF_temp_autumn[1] ~ dnorm(0, 1)
  beta_PTF_grass_autumn[1] <- 0
  beta_PTF_ag_autumn[1] <- 0
  beta_PTF_bog_autumn[1] ~ dnorm(0, 1)
  
  ## effects of behavior on survival
  # random intercept for individual
  for(i in ind_autumn[]){ # individual index
    alpha_surv_autumn[i] ~ dnorm(0, tau_surv_autumn)
  } # i
  # individual-specific variation
  tau_surv_autumn <- pow(sd_surv_autumn, -2)
  sd_surv_autumn ~ dunif(0, 10)
  
  # mean betas
  mu_beta_ODBA_autumn ~ dnorm(0, 0.01)
  mu_beta_PTF_autumn ~ dnorm(0, 0.01)
  
  # sub_seasonal betas
  for(s in c(1, 2, 3)){ # breeding, first migratory flight, first half staging
    beta_ODBA_autumn[s] ~ dnorm(mu_beta_ODBA_autumn, tau_beta_ODBA_autumn)
  } # s
  tau_beta_ODBA_autumn <- pow(sd_beta_ODBA_autumn, -2)
  sd_beta_ODBA_autumn ~ dunif(0, 10)
  
  for(s in c(1, 3)){ # breeding, first half staging
    beta_PTF_autumn[s] ~ dnorm(mu_beta_PTF_autumn, tau_beta_PTF_autumn)
  } # s
  tau_beta_PTF_autumn <- pow(sd_beta_PTF_autumn, -2)
  sd_beta_PTF_autumn ~ dunif(0, 10)
  
  # departure date beta
  beta_depart_breed ~ dnorm(0, 0.01)
  
  # random year effect
  for(t in 1:nyears){ # year index
    eps_year_autumn[t] ~ dnorm(0, tau_year_autumn)
  } # t
  tau_year_autumn <- pow(sd_year_autumn, -2)
  sd_year_autumn ~ dunif(0, 10)
  
  #### likelihood
  ### spring model
  ## effects of weather and habitat use on behavior
  # observation variation
  for(b in 1:nobs2){ # individual/year/sub-season index
    log_ODBA_spring[b] ~ dnorm(mu_log_ODBA_spring[ids2[b], years2[b], sub_seasons2[b]], 
                               tau_ODBA_obs_spring[ids2[b],years2[b]])
  
  # link function
  mu_log_ODBA_spring[ids2[b], years2[b], sub_seasons2[b]] <- 
      alpha_ODBA_spring[ids2[b]] + 
      beta_ODBA_precip_spring[sub_seasons2[b]] * precip_spring2[b] + 
      beta_ODBA_temp_spring[sub_seasons2[b]] * temp_spring2[b] + 
      beta_ODBA_storm_spring[sub_seasons2[b]] * storm_spring2[b] +
      beta_ODBA_grass_spring[sub_seasons2[b]] * grass_spring2[b] + 
      beta_ODBA_ag_spring[sub_seasons2[b]] * ag_spring2[b] + 
      beta_ODBA_bog_spring[sub_seasons2[b]] * bog_spring2[b]
  } # b
  
  # observation variation
  for(c in 1:nobs3){ # individual/year/(non-migratory) sub-season index
    PTF_spring[c] ~ dbin(p_PTF_spring[ids3[c], years3[c], sub_seasons3[c]], n_spring[c])
      
  # link function
    logit(p_PTF_spring[ids3[c], years3[c], sub_seasons3[c]]) <- 
      alpha_PTF_spring[ids3[c]] + 
      beta_PTF_precip_spring[sub_seasons3[c]] * precip_spring3[c] + 
      beta_PTF_temp_spring[sub_seasons3[c]] * temp_spring3[c] +
      beta_PTF_grass_spring[sub_seasons3[c]] * grass_spring3[c] + 
      beta_PTF_ag_spring[sub_seasons3[c]] * ag_spring3[c] + 
      beta_PTF_bog_spring[sub_seasons3[c]] * bog_spring3[c]
  } # c
  
  ## effects of behavior on breeding outcome
  for(a in 1:nobs1){ # individual/year index
    B_spring[a] ~ dcat(p[a, 1:n_B])
    
    for(k in 1:n_B){ # breeding outcome index
  # link function
      logit(q[a, k]) <- 
        alpha_B_spring[ids1[a], k] + 
        beta_ODBA_spring[1, k] * mu_log_ODBA_spring[ids1[a], years1[a], 1] +
        beta_ODBA_spring[2, k] * mu_log_ODBA_spring[ids1[a], years1[a], 2] +
        beta_ODBA_spring[3, k] * mu_log_ODBA_spring[ids1[a], years1[a], 3] +
        beta_ODBA_spring[4, k] * mu_log_ODBA_spring[ids1[a], years1[a], 4] +
        beta_ODBA_spring[5, k] * mu_log_ODBA_spring[ids1[a], years1[a], 5] +
        beta_ODBA_spring[6, k] * mu_log_ODBA_spring[ids1[a], years1[a], 6] +
        beta_PTF_spring[1, k] *  logit(p_PTF_spring[ids1[a], years1[a], 1]) +
        beta_PTF_spring[3, k] *  logit(p_PTF_spring[ids1[a], years1[a], 3]) +
        beta_PTF_spring[4, k] *  logit(p_PTF_spring[ids1[a], years1[a], 4]) +
        beta_PTF_spring[6, k] *  logit(p_PTF_spring[ids1[a], years1[a], 6]) +
        beta_arrival_breed[k] * arrival_day_breed[a] +
        eps_year_spring[years1[a], k]
  # transform into probabilities that sum to 1
      p[a, k] <- q[a, k] / sum(q[a, 1:n_B])
    } # k
  } # a
  
  ### autumn model
  ## effects of breeding outcome, weather, and habitat use on behavior
  # observation variation
  for(e in 1:nobs5){ # individual/year/sub-season index
    log_ODBA_autumn[e] ~ dnorm(mu_log_ODBA_autumn[ids5[e], years5[e], sub_seasons5[e]], 
                               tau_ODBA_obs_autumn[ids5[e],years5[e]])
  
  # link function
    mu_log_ODBA_autumn[ids5[e], years5[e], sub_seasons5[e]] <- 
      alpha_ODBA_autumn[ids5[e]] + 
      beta_ODBA_B_autumn[sub_seasons5[e], B_autumn5[e]] + 
      beta_ODBA_precip_autumn[sub_seasons5[e]] * precip_autumn5[e] + 
      beta_ODBA_temp_autumn[sub_seasons5[e]] * temp_autumn5[e] + 
      beta_ODBA_storm_autumn[sub_seasons5[e]] * storm_autumn5[e] +
      beta_ODBA_grass_autumn[sub_seasons5[e]] * grass_autumn5[e] + 
      beta_ODBA_ag_autumn[sub_seasons5[e]] * ag_autumn5[e] + 
      beta_ODBA_bog_autumn[sub_seasons5[e]] * bog_autumn5[e]
  } # e
  
  # observation variation
  for(f in 1:nobs6){ # individual/year/(non-migratory) sub-season index
    PTF_autumn[f] ~ dbin(p_PTF_autumn[ids6[f], years6[f], sub_seasons6[f]], n_autumn[f])
    
    # link function
    logit(p_PTF_autumn[ids6[f], years6[f], sub_seasons6[f]]) <- 
      alpha_PTF_autumn[ids6[f]] + 
      beta_PTF_B_autumn[sub_seasons6[f], B_autumn6[f]] + 
      beta_PTF_precip_autumn[sub_seasons6[f]] * precip_autumn6[f] + 
      beta_PTF_temp_autumn[sub_seasons6[f]] * temp_autumn6[f] +
      beta_PTF_grass_autumn[sub_seasons6[f]] * grass_autumn6[f] + 
      beta_PTF_ag_autumn[sub_seasons6[f]] * ag_autumn6[f] +
      beta_PTF_bog_autumn[sub_seasons6[f]] * bog_autumn6[f]
  } # f
  
  ## effects of behavior on survival
  for(d in 1:nobs4){ # individual/year index
    surv_autumn[d] ~ dbern(phi[ids4[d], years4[d]])
    
    # link function
    logit(phi[ids4[d], years4[d]]) <- alpha_surv_autumn[ids4[d]] + 
      beta_ODBA_autumn[1] * mu_log_ODBA_autumn[ids4[d], years4[d], 1] + 
      beta_ODBA_autumn[2] * mu_log_ODBA_autumn[ids4[d], years4[d], 2] + 
      beta_ODBA_autumn[3] * mu_log_ODBA_autumn[ids4[d], years4[d], 3] + 
      beta_PTF_autumn[1] * logit(p_PTF_autumn[ids4[d], years4[d], 1]) +
      beta_PTF_autumn[3] * logit(p_PTF_autumn[ids4[d], years4[d], 3]) +
      beta_depart_breed * depart_day_breed[d] + 
      eps_year_autumn[years4[d]]
  } # d
}")

# set data
jags_data <- list(
  # i loop (full model, individual index)
  nind = length(unique(c(spring_dat_all$id, autumn_dat_all$id))),
  ind_autumn = sort(unique(autumn_dat_all$id)),
  # k loop (full model, breeding outcome index)
  n_B = length(unique(spring_dat_all$breeding_outcome)),
  # s loop (full model, sub-season index)
  n_sub_seasons = length(unique(spring_dat_all$sub_season)),
  # t loop (full model, year index)
  nyears = length(unique(spring_dat_all$year)),
  # 1. a loop (spring model, individual/year index)
  nobs1 = nrow(spring_dat_annual), ids1 = spring_dat_annual$id, 
  years1 = spring_dat_annual$year,
  B_spring = spring_dat_annual$breeding_outcome,
  arrival_day_breed = spring_dat_annual$scaled_arrival_day,
  # 2. b loop (spring model, individual/year/sub-season index)
  nobs2 = nrow(spring_dat_all), ids2 = spring_dat_all$id, 
  years2 = spring_dat_all$year, sub_seasons2 = spring_dat_all$sub_season,
  log_ODBA_spring = spring_dat_all$log_ODBA,
  precip_spring2 = spring_dat_all$scaled_mean_precip, 
  temp_spring2 = spring_dat_all$scaled_prop_days_below_freezing, 
  storm_spring2 = spring_dat_all$scaled_prop_storm_days,
  grass_spring2 = spring_dat_all$scaled_prop_grass,
  ag_spring2 = spring_dat_all$scaled_prop_ag,
  bog_spring2 = spring_dat_all$scaled_prop_bog,
  # 3. c loop (spring model, individual/year/(non-migratory) sub-season index)
  n_spring = spring_dat_non_mig$num_ACC_fixes,
  nobs3 = nrow(spring_dat_non_mig), ids3 = spring_dat_non_mig$id, 
  years3 = spring_dat_non_mig$year, sub_seasons3 = spring_dat_non_mig$sub_season,
  precip_spring3 = spring_dat_non_mig$scaled_mean_precip, 
  temp_spring3 = spring_dat_non_mig$scaled_prop_days_below_freezing, 
  grass_spring3 = spring_dat_non_mig$scaled_prop_grass,
  ag_spring3 = spring_dat_non_mig$scaled_prop_ag,
  bog_spring3 = spring_dat_non_mig$scaled_prop_bog,
  PTF_spring = spring_dat_non_mig$num_feed_fixes,
  # 4. d loop (autumn model, individual/year index)
  nobs4 = nrow(autumn_dat_annual), ids4 = autumn_dat_annual$id, 
  years4 = autumn_dat_annual$year, surv_autumn = autumn_dat_annual$annual_survival,
  depart_day_breed = autumn_dat_annual$scaled_departure_day,
  # 5. e loop (autumn model, individual/year/sub-season index)
  nobs5 = nrow(autumn_dat_all), ids5 = autumn_dat_all$id, 
  years5 = autumn_dat_all$year, sub_seasons5 = autumn_dat_all$sub_season,
  log_ODBA_autumn = autumn_dat_all$log_ODBA,
  precip_autumn5 = autumn_dat_all$scaled_mean_precip, 
  temp_autumn5 = autumn_dat_all$scaled_prop_days_below_freezing, 
  storm_autumn5 = autumn_dat_all$scaled_prop_storm_days,
  grass_autumn5 = autumn_dat_all$scaled_prop_grass,
  ag_autumn5 = autumn_dat_all$scaled_prop_ag,
  bog_autumn5 = autumn_dat_all$scaled_prop_bog,
  B_autumn5 = autumn_dat_all$breeding_outcome,
  # 6. f loop (autumn model, individual/year/(non-migratory) sub-season index)
  n_autumn = autumn_dat_non_mig$num_ACC_fixes,
  nobs6 = nrow(autumn_dat_non_mig), ids6 = autumn_dat_non_mig$id, 
  years6 = autumn_dat_non_mig$year, sub_seasons6 = autumn_dat_non_mig$sub_season,
  precip_autumn6 = autumn_dat_non_mig$scaled_mean_precip, 
  temp_autumn6 = autumn_dat_non_mig$scaled_prop_days_below_freezing, 
  grass_autumn6 = autumn_dat_non_mig$scaled_prop_grass,
  ag_autumn6 = autumn_dat_non_mig$scaled_prop_ag,
  bog_autumn6 = autumn_dat_non_mig$scaled_prop_bog,
  PTF_autumn = autumn_dat_non_mig$num_feed_fixes,
  B_autumn6 = autumn_dat_non_mig$breeding_outcome
)

# set initial values
annual_inits <- function(dat){
  a_inits <- unique(dat %>% select(id, year))
  a_inits$init <- runif(nrow(a_inits), 0, 10)
  reshape2::acast(a_inits, id ~ year, value.var = "init", na.rm = F)
}

annual_inits2 <- function(dat1, dat2, nyears){
  a_inits <- unique(dat1 %>% select(id, year))
  a_inits$init <- runif(nrow(a_inits), 0, 10)
  a_inits2 <- data.frame(id = rep(1:max(dat2$id), times = nyears), 
                         year = rep(1:nyears, each = max(dat2$id)))
  a_inits <- merge(a_inits2, a_inits, by = c("id", "year"), all.x = T)
  a_inits <- reshape2::acast(a_inits, id ~ year, value.var = "init", na.rm = F, drop = F)
  a_inits
}

na_fill_inits <- function(n_param, na_param, sd){
  t_inits <- rnorm(n_param, mean = 0, sd = sd)
  t_inits[na_param] <- NA
  t_inits
}

na_fill_inits2 <- function(n_sub_seasons, n_B, na_sub_seasons, na_B, sd){
  t_inits <- matrix(rnorm(n_sub_seasons * n_B, mean = 0, sd = sd), nrow = n_sub_seasons, ncol = n_B)
  if(any(is.na(na_sub_seasons) == F)){
    t_inits[na_sub_seasons, ] <- NA
  }
  t_inits[ , na_B] <- NA
  t_inits
}

na_fill_inits3 <- function(n_param, na_param, max){
  t_inits <- runif(n_param, 0, max)
  t_inits[na_param] <- NA
  t_inits
}

inits <- function() {
  list(
    sd_ODBA_spring = runif(1, 0, 1),
    sd_PTF_spring = runif(1, 0, 1),
    sd_ODBA_obs_spring = annual_inits(spring_dat_all),
    beta_ODBA_precip_spring = rnorm(6, mean = 0, sd = 1),
    beta_ODBA_temp_spring = rnorm(6, mean = 0, sd = 1),
    beta_ODBA_storm_spring = na_fill_inits(6, c(1, 3, 4, 6), 1),
    beta_ODBA_grass_spring = na_fill_inits(6, c(2, 5, 6), 1),
    beta_ODBA_ag_spring = na_fill_inits(6, c(2, 5, 6), 1),
    beta_ODBA_bog_spring = na_fill_inits(6, c(2, 5), 1),
    beta_PTF_precip_spring = na_fill_inits(6, c(2, 5), 1),
    beta_PTF_temp_spring = na_fill_inits(6, c(2, 5), 1),
    beta_PTF_grass_spring = na_fill_inits(6, c(2, 5, 6), 1),
    beta_PTF_ag_spring = na_fill_inits(6, c(2, 5, 6), 1),
    beta_PTF_bog_spring = na_fill_inits(6, c(2, 5), 1),
    sd_B_spring = runif(2, 1, 10),
    mu_beta_ODBA_spring = rnorm(2, 0, 1),
    sd_beta_ODBA_spring = runif(2, 1, 10),
    mu_beta_PTF_spring = rnorm(2, 0, 1),
    sd_beta_PTF_spring = runif(2, 1, 10),
    beta_arrival_breed = na_fill_inits(3, 3, 1),
    sd_year_spring = runif(2, 1, 10),
    sd_ODBA_autumn = runif(1, 0, 1),
    sd_PTF_autumn = runif(1, 0, 1),
    sd_ODBA_obs_autumn = annual_inits2(autumn_dat_all, spring_dat_all, nyears = 5)[-49,],
    beta_ODBA_B_autumn = na_fill_inits2(6, 3, NA, 1, 1),
    beta_ODBA_precip_autumn = rnorm(6, mean = 0, sd = 1),
    beta_ODBA_temp_autumn = rnorm(6, mean = 0, sd = 1),
    beta_ODBA_storm_autumn = na_fill_inits(6, c(1, 3, 4, 6), 1),
    beta_ODBA_grass_autumn = na_fill_inits(6, c(1, 2, 5), 1),
    beta_ODBA_ag_autumn = na_fill_inits(6, c(1, 2, 5), 1),
    beta_ODBA_bog_autumn = na_fill_inits(6, c(2, 5), 1),
    beta_PTF_B_autumn = na_fill_inits2(6, 3, c(2, 5), 1, 1),
    beta_PTF_precip_autumn = na_fill_inits(6, c(2, 5), 1),
    beta_PTF_temp_autumn = na_fill_inits(6, c(2, 5), 1),
    beta_PTF_grass_autumn = na_fill_inits(6, c(1, 2, 5), 1),
    beta_PTF_ag_autumn = na_fill_inits(6, c(1, 2, 5), 1),
    beta_PTF_bog_autumn = na_fill_inits(6, c(2, 5), 1),
    sd_surv_autumn = runif(1, 0, 10),
    mu_beta_ODBA_autumn = rnorm(1, 0, 1),
    sd_beta_ODBA_autumn = runif(1, 0, 10),
    mu_beta_PTF_autumn = rnorm(1, 0, 1),
    sd_beta_PTF_autumn = runif(1, 0, 10),
    beta_depart_breed = rnorm(1, 0, 1),
    sd_year_autumn = runif(1, 0, 10)
  )
}

params <- c("alpha_ODBA_spring", 
            "sd_ODBA_spring",
            "alpha_PTF_spring",
            "sd_PTF_spring",
            "sd_ODBA_obs_spring",
            "beta_ODBA_precip_spring", 
            "beta_ODBA_temp_spring", 
            "beta_ODBA_storm_spring",
            "beta_ODBA_grass_spring",
            "beta_ODBA_ag_spring", 
            "beta_ODBA_bog_spring",
            "beta_PTF_precip_spring",
            "beta_PTF_temp_spring",
            "beta_PTF_grass_spring",
            "beta_PTF_ag_spring",
            "beta_PTF_bog_spring",
            "alpha_B_spring",
            "sd_B_spring",
            "beta_ODBA_spring",
            "beta_PTF_spring",
            "mu_beta_ODBA_spring",
            "sd_beta_ODBA_spring",
            "mu_beta_PTF_spring",
            "sd_beta_PTF_spring",
            "beta_arrival_breed",
            "sd_year_spring",
            "alpha_ODBA_autumn",
            "sd_ODBA_autumn",
            "alpha_PTF_autumn",
            "sd_PTF_autumn",
            "sd_ODBA_obs_autumn",
            "beta_ODBA_B_autumn",
            "beta_ODBA_precip_autumn",
            "beta_ODBA_temp_autumn",
            "beta_ODBA_storm_autumn",
            "beta_ODBA_grass_autumn",
            "beta_ODBA_ag_autumn",
            "beta_ODBA_bog_autumn",
            "beta_PTF_B_autumn",
            "beta_PTF_precip_autumn",
            "beta_PTF_temp_autumn",
            "beta_PTF_grass_autumn",
            "beta_PTF_ag_autumn",
            "beta_PTF_bog_autumn",
            "alpha_surv_autumn",
            "sd_surv_autumn",
            "beta_ODBA_autumn",
            "beta_PTF_autumn",
            "mu_beta_ODBA_autumn",
            "sd_beta_ODBA_autumn",
            "mu_beta_PTF_autumn",
            "sd_beta_PTF_autumn",
            "beta_depart_breed",
            "sd_year_autumn")

# Run model
library(jagsUI)
out <- autojags(jags_data, inits, params, "full_annual_cycle_model_categorical.txt", 
                n.chains = 3, iter.increment = 12000, n.burnin = 20000, n.thin = 2, Rhat.limit = 1.1, max.iter = 300000)

