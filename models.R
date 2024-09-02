#Total biomass model----
T_mini_model <- lmer(log_totbio ~ log_ht*myc_group + (1|study_species), data = full_df_mod)
summary(T_mini_model)

#test full model with all interaction terms
T_full_model <-lmer(log_totbio ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(T_full_model)
summary(T_full_model)

library(brms)
model <- brm(bf(log(m.to) ~ log_ht*myc_group + (1|study_species)), 
             data = full_df_mod, 
             family = gaussian(), 
             chains = 4, 
             iter = 4000, 
             warmup = 2000, 
             cores = 4, 
             control = list(adapt_delta = 0.9),
             backend = "cmdstanr",
             seed = 123)
summary(model)
pp_check(model, ndraws = 100)
conditional_effects(model)




#next, test full model with all interaction terms
R_full_model <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(R_full_model)
summary(R_full_model)

full_df_mod$RSL = cbind(roots = full_df_mod$roots,
                        shoots = full_df_mod$shoots,
                        leaves = full_df_mod$leaves)

R_full_model <- brm(bf(RSL ~ log_ht*myc_group), 
             data = full_df_mod, 
             family = dirichlet(), 
             chains = 4, 
             iter = 6000, 
             warmup = 3000, 
             cores = 4, 
             control = list(adapt_delta = 0.95),
             backend = "cmdstanr",
             seed = 123)
summary(R_full_model)
pp_check(R_full_model, ndraws = 100)
conditional_effects(R_full_model,
                    categorical = TRUE)
