within.sp = full_df_mod %>% 
  rename(roots = RmTm,
         stems = SmTm,
         leaves = LmTm) %>% 
  mutate(RSL = cbind(roots = .$roots,
                     stems = .$stems,
                     leaves = .$leaves)) %>% 
  group_by(SppName) %>% 
  mutate(.after = h.t,
         ht.median = median(h.t),
         rel.ht = log2(h.t/ht.median),
         N = n()) %>% 
  ungroup()



library(brms)
library(emmeans)

m.Bet_alle.s <- brm(bf(stems ~ rel.ht + (1|Study)), 
                    data = within.sp %>% 
                      filter(SppName == "Betula_alleghaniensis"), 
                    family = Beta(), 
                    chains = 4, 
                    iter = 6000, 
                    warmup = 3000, 
                    cores = 4, 
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 12),
                    backend = "cmdstanr",
                    file = "fits/m.Bet_alle.s",
                    seed = 123)
pp_check(m.Bet_alle.s, ndraws = 200)
plot(conditional_effects(m.Bet_alle.s),
     points = TRUE)

m.Bet_alle.r <- brm(bf(roots ~ rel.ht + (1|Study)), 
                    data = within.sp %>% 
                      filter(SppName == "Betula_alleghaniensis"), 
                    family = Beta(), 
                    chains = 4, 
                    iter = 6000, 
                    warmup = 3000, 
                    cores = 4, 
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    file = "fits/m.Bet_alle.r",
                    seed = 123)
pp_check(m.Bet_alle.r, ndraws = 200)
plot(conditional_effects(m.Bet_alle.r),
     points = TRUE)


m.Bet_alle.l <- brm(bf(leaves ~ rel.ht + (1|Study)), 
                    data = within.sp %>% 
                      filter(SppName == "Betula_alleghaniensis"), 
                    family = Beta(), 
                    chains = 4, 
                    iter = 6000, 
                    warmup = 3000, 
                    cores = 4, 
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 12),
                    backend = "cmdstanr",
                    file = "fits/m.Bet_alle.l",
                    seed = 123)
pp_check(m.Bet_alle.l, ndraws = 200)
plot(conditional_effects(m.Bet_alle.l),
     points = TRUE)



m.Bet_alle.rsl <- brm(bf(RSL ~ rel.ht + (1|Study)), 
                      data = within.sp %>% 
                        filter(SppName == "Betula_alleghaniensis"), 
                      family = dirichlet(), 
                      chains = 4, 
                      iter = 9000, 
                      warmup = 3000, 
                      cores = 4, 
                      control = list(adapt_delta = 0.99),
                      backend = "cmdstanr",
                      file = "fits/m.Bet_alle.rsl",
                      seed = 123)
summary(m.Bet_alle.rsl)
pp_check(m.Bet_alle.rsl, ndraws = 100)
conditional_effects(m.Bet_alle.rsl,
                    categorical = TRUE)







# species with at least 10 individuals
m.most.s <- brm(bf(stems ~ rel.ht + (1|SppName/Study)), 
                data = within.sp %>% 
                  filter(N >= 10), 
                family = Beta(), 
                chains = 4, 
                iter = 6000, 
                warmup = 3000, 
                cores = 4, 
                control = list(adapt_delta = 0.95),
                backend = "cmdstanr",
                file = "fits/m.most.s",
                seed = 123)
