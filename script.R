############################ Roots,stems & leaves ##############################

# author: Angelos Amyntas
# date: 20240907

# This script contains the complete code for the analysis presented in 
# [https://amynang.github.io/posts/20240907-roots-stems-leaves/]

########################### load data, standardize height ######################
library(tidyverse)

full_df_mod_new = read.csv("data/full_df_mod_new.csv")

within.sp = full_df_mod_new %>% 
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


######################## Beta vs Dirichlet for one species #####################
library(brms)

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
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 12),
                      backend = "cmdstanr",
                      file = "fits/m.Bet_alle.rsl",
                      seed = 123)
summary(m.Bet_alle.rsl)
pp_check(m.Bet_alle.rsl, ndraws = 200)
conditional_effects(m.Bet_alle.rsl,
                    categorical = TRUE)

#################### posterior predictive check for Dirichlet ##################
library(bayesplot)
# generate predictions
yrep <- posterior_predict(m.Bet_alle.rsl, ndraws = 200)
# yrep is a three dimensional array:
# 200 predictions x 295 observations x 3 response components
dim(yrep)
# here's one slice:
yrep[1, 1, 1:3]
#     roots     stems    leaves 
# 0.3791285 0.4160540 0.2048174

# three custom-made density overlay plots, similar to the default pp_check() output
ppc_dens_overlay(within.sp$roots[within.sp$SppName == "Betula_alleghaniensis"], 
                 yrep[ , 1:295, 1])
ppc_dens_overlay(within.sp$stems[within.sp$SppName == "Betula_alleghaniensis"], 
                 yrep[ , 1:295, 2])
ppc_dens_overlay(within.sp$leaves[within.sp$SppName == "Betula_alleghaniensis"], 
                 yrep[ , 1:295, 3])


################### Plot model estimates for one species #######################
library(tidybayes)
library(modelr)

root = within.sp %>%
  data_grid(rel.ht = seq_range(rel.ht, n = 51)) %>%
  add_epred_draws(m.Bet_alle.r, re_formula = NA) %>% 
  add_column(.category = "roots")

stem = within.sp %>%
  data_grid(rel.ht = seq_range(rel.ht, n = 51)) %>%
  add_epred_draws(m.Bet_alle.s, re_formula = NA) %>% 
  add_column(.category = "stems")

leaf = within.sp %>%
  data_grid(rel.ht = seq_range(rel.ht, n = 51)) %>%
  add_epred_draws(m.Bet_alle.l, re_formula = NA) %>% 
  add_column(.category = "leaves")

rsl = rbind(root,stem,leaf)
rsl$.category = factor(rsl$.category, levels = c("roots","stems","leaves"))

fig1a = rsl %>%
  ggplot(aes(x = rel.ht, y = .epred, color = .category,
             fill = .category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_lineribbon(aes(y = .epred), .width = c(0.9),
                  point_interval = "median_qi",
                  linewidth = .2)+
  scale_x_continuous(breaks = c(-3,-1,0,1,3), 
                     labels = c("1/8","1/2",1,2,8)) +
  scale_color_manual(values = c("#efd08e","#8d8067","#9dbc9b")) +
  scale_fill_manual(values = c("#efd08e80","#8d806780","#9dbc9b80")) +
  ylab("Biomass proportion") +
  xlab("multiplicative change in relative height") +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Predictions of three Beta regressions")

fig1b = within.sp %>%
  data_grid(rel.ht = seq_range(rel.ht, n = 51)) %>%
  add_epred_draws(m.Bet_alle.rsl, re_formula = NA) %>%
  ggplot(aes(x = rel.ht, y = .epred, color = .category,
             fill = .category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_lineribbon(aes(y = .epred), .width = c(0.9),
                  point_interval = "median_qi",
                  linewidth = .2)+
  scale_x_continuous(breaks = c(-3,-1,0,1,3), 
                     labels = c("1/8","1/2",1,2,8)) +
  scale_color_manual(values = c("#efd08e","#8d8067","#9dbc9b")) +
  scale_fill_manual(values = c("#efd08e80","#8d806780","#9dbc9b80")) +
  ylab("Biomass proportion") +
  xlab("multiplicative change in relative height") +
  theme_bw() +
  theme(legend.title = element_blank()) + 
  ggtitle("Predictions of Dirichlet regression")

library(patchwork)
fig1 = fig1a+fig1b

ggsave("figs/fig1.png", fig1, units = "mm",
       width = 180,
       height = 75,
       scale = 1.4)



############################ three Betas for all species #######################

m.all.r_hhm2 <- brm(bf(roots ~ (rel.ht  + leaf_habit + myc_group)^2 + 
                         (1 + rel.ht|SppName) + (1|Study)), 
                    data = within.sp, 
                    family = Beta(), 
                    chains = 4, 
                    iter = 9000, 
                    warmup = 3000, 
                    cores = 4, 
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 12),
                    backend = "cmdstanr",
                    file = "fits/m.all.r_hhm2",
                    seed = 123)
summary(m.all.r_hhm2)
pp_check(m.all.r_hhm2, ndraws = 200)

m.all.s_hhm2 <- brm(bf(stems ~ (rel.ht  + leaf_habit + myc_group)^2 + 
                         (1 + rel.ht|SppName) + (1|Study)), 
                    data = within.sp, 
                    family = Beta(), 
                    chains = 4, 
                    iter = 9000, 
                    warmup = 3000, 
                    cores = 4, 
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 12),
                    backend = "cmdstanr",
                    file = "fits/m.all.s_hhm2",
                    seed = 123)
summary(m.all.s_hhm2)
pp_check(m.all.s_hhm2, ndraws = 200)

m.all.l_hhm2 <- brm(bf(leaves ~ (rel.ht  + leaf_habit + myc_group)^2 + 
                         (1 + rel.ht|SppName) + (1|Study)), 
                    data = within.sp, 
                    family = Beta(), 
                    chains = 4, 
                    iter = 9000, 
                    warmup = 3000, 
                    cores = 4, 
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 12),
                    backend = "cmdstanr",
                    file = "fits/m.all.l_hhm2",
                    seed = 123)
summary(m.all.l_hhm2)
pp_check(m.all.l_hhm2, ndraws = 200)


############################ Dirichlet for all species #########################

m.all.rsl_hhm2 <- brm(bf(RSL ~ (rel.ht  + leaf_habit + myc_group)^2 + 
                           (1 + rel.ht|SppName) + (1|Study)), 
                      data = within.sp, 
                      family = dirichlet(), 
                      chains = 4, 
                      iter = 9000, 
                      warmup = 3000, 
                      cores = 4, 
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 12),
                      backend = "cmdstanr",
                      file = "fits/m.all.rsl_hhm2",
                      seed = 123)
summary(m.all.rsl_hhm2)
# generate predictions
yrep <- posterior_predict(m.all.rsl_hhm2, ndraws = 200)
# three custom-made density overlay plots, similar to the default pp_check() output
ppc_dens_overlay(within.sp$roots, 
                 yrep[ , 1:1429, 1])
ppc_dens_overlay(within.sp$stems, 
                 yrep[ , 1:1429, 2])
ppc_dens_overlay(within.sp$leaves, 
                 yrep[ , 1:1429, 3])


fig2 = within.sp %>%
  data_grid(rel.ht = seq_range(rel.ht, n = 51),
            leaf_habit = c("evergreen", 
                           "deciduous"),
            myc_group = c("ECM","AM")) %>%
  add_epred_draws(m.all.rsl_hhm2, re_formula = NA) %>%
  ggplot(aes(x = rel.ht, y = .epred, color = .category,
             fill = .category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_lineribbon(aes(y = .epred), .width = c(0.9),
                  point_interval = "median_qi",
                  linewidth = .2)+
  scale_x_continuous(breaks = c(-3,-1,0,1,3), 
                     labels = c("1/8","1/2",1,2,8)) +
  facet_grid(myc_group~leaf_habit) + 
  scale_color_manual(values = c("#efd08e","#8d8067","#9dbc9b")) +
  scale_fill_manual(values = c("#efd08e80","#8d806780","#9dbc9b80")) +
  ylab("Biomass proportion") +
  xlab("multiplicative change in relative height") +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("figs/fig2.png", fig2, units = "mm",
       width = 180,
       height = 140,
       scale = 1.4)


fig3a = within.sp %>% 
  data_grid(expand_grid(rel.ht = c(0),
                        leaf_habit = c("deciduous"),
                        myc_group = c("ECM","AM"))) %>% 
  add_epred_draws(m.all.rsl_hhm2,
                  re_formula = NA) %>% 
  mutate(group = str_c(.category, myc_group)) %>% 
  ggplot(aes(x = .epred, fill = factor(group))) +
  stat_halfeye(aes(shape = myc_group),
               .width = c(.9),
               size = 5,
               alpha = .5) +
  scale_fill_manual(values = c("#9dbc9b","#9dbc9b",
                               "#efd08e","#efd08e",
                               "#8d8067","#8d8067"),
                    guide = "none") +
  ylab("") +
  xlab("estimate") +
  facet_wrap(~.category,
             scales = "free") +
  theme_bw() + 
  ggtitle("Biomass allocation of deciduous trees at median height")

fig3b = within.sp %>% 
  data_grid(expand_grid(rel.ht = c(0),
                        leaf_habit = c("evergreen"),
                        myc_group = c("ECM","AM"))) %>% 
  add_epred_draws(m.all.rsl_hhm2,
                  re_formula = NA) %>% 
  mutate(group = str_c(.category, myc_group)) %>% 
  ggplot(aes(x = .epred, fill = factor(group))) +
  stat_halfeye(aes(shape = myc_group),
               .width = c(.9),
               size = 5,
               alpha = .5) +
  scale_fill_manual(values = c("#9dbc9b","#9dbc9b",
                               "#efd08e","#efd08e",
                               "#8d8067","#8d8067"),
                    guide = "none") +
  ylab("") +
  xlab("estimate") +
  facet_wrap(~.category,
             scales = "free") +
  theme_bw() + 
  ggtitle("Biomass allocation of evergreen trees at median height")


fig3 = fig3a/fig3b
ggsave("figs/fig3.png", fig3, units = "mm",
       width = 180,
       height = 120,
       scale = 1.4)


############################# Phylogenetic regression ##########################

# see relevant brms vignette
# https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html#a-phylogenetic-model-with-repeated-measurements

library(ape)
# load the tree
phylo = read.tree("data/phyliptree1.phy")
# phylogenetic variance-covariance matrix
A = ape::vcv.phylo(phylo)
# replaces "name" with name
dimnames(A) = list(unique(within.sp$SppName),
                   unique(within.sp$SppName))
# grouping of observations by species (brms requires a distinct named variable)
within.sp$phylo = within.sp$SppName

m.all.rsl_hhm2_h.ph <- brm(bf(RSL ~ (rel.ht  + leaf_habit + myc_group)^2 + 
                                (1 + rel.ht|gr(phylo, cov = A)) + (1 + rel.ht|SppName) + (1|Study)), 
                           data = within.sp, 
                           data2 = list(A = A),
                           family = dirichlet(), 
                           chains = 4, 
                           iter = 9000, 
                           warmup = 3000, 
                           cores = 4, 
                           control = list(adapt_delta = 0.99,
                                          max_treedepth = 12),
                           backend = "cmdstanr",
                           file = "fits/m.all.rsl_hhm2_h.ph",
                           seed = 123)

summary(m.all.rsl_hhm2_h.ph)
# generate predictions
yrep <- posterior_predict(m.all.rsl_hhm2_h.ph, ndraws = 200)
# three custom-made density overlay plots, similar to the default pp_check() output
ppc_dens_overlay(within.sp$roots, 
                 yrep[ , 1:1429, 1])
ppc_dens_overlay(within.sp$stems, 
                 yrep[ , 1:1429, 2])
ppc_dens_overlay(within.sp$leaves, 
                 yrep[ , 1:1429, 3])


fig4 = within.sp %>%
  data_grid(rel.ht = seq_range(rel.ht, n = 51),
            leaf_habit = c("evergreen", 
                           "deciduous"),
            myc_group = c("ECM","AM")) %>%
  add_epred_draws(m.all.rsl_hhm2_h.ph, re_formula = NA) %>%
  ggplot(aes(x = rel.ht, y = .epred, color = .category,
             fill = .category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_lineribbon(aes(y = .epred), .width = c(0.9),
                  point_interval = "median_qi",
                  linewidth = .2)+
  scale_x_continuous(breaks = c(-3,-1,0,1,3), 
                     labels = c("1/8","1/2",1,2,8)) +
  facet_grid(myc_group~leaf_habit) + 
  scale_color_manual(values = c("#efd08e","#8d8067","#9dbc9b")) +
  scale_fill_manual(values = c("#efd08e80","#8d806780","#9dbc9b80")) +
  ylab("Biomass proportion") +
  xlab("multiplicative change in relative height") +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("figs/fig4.png", fig4, units = "mm",
       width = 180,
       height = 140,
       scale = 1.4)












m.ac.all.rsl_hhm2_h.ph <- brm(bf(RSL ~ (rel.ht2  + leaf_habit + myc_group)^2 + 
                                   (1 + rel.ht2|gr(phylo, cov = A)) + (1 + rel.ht2|SppName) + (1|Study)), 
                              data = across.sp, 
                              data2 = list(A = A),
                              family = dirichlet(), 
                              chains = 4, 
                              iter = 9000, 
                              warmup = 3000, 
                              cores = 4, 
                              control = list(adapt_delta = 0.99,
                                             max_treedepth = 12),
                              backend = "cmdstanr",
                              file = "fits/m.ac.all.rsl_hhm2_h.ph",
                              seed = 123)

within.sp %>%
  data_grid(rel.ht2 = seq_range(rel.ht, n = 51),
            leaf_habit = c("evergreen", 
                           "deciduous"),
            myc_group = c("ECM","AM")) %>%
  add_epred_draws(m.ac.all.rsl_hhm2_h.ph, re_formula = NA) %>%
  ggplot(aes(x = rel.ht2, y = .epred, color = .category,
             fill = .category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_lineribbon(aes(y = .epred), .width = c(0.9),
                  point_interval = "median_qi",
                  linewidth = .2)+
  #scale_x_continuous(breaks = c(-3,-1,0,1,3), 
  #                   labels = c("1/8","1/2",1,2,8)) +
  facet_grid(myc_group~leaf_habit) + 
  scale_color_manual(values = c("#efd08e","#8d8067","#9dbc9b")) +
  scale_fill_manual(values = c("#efd08e80","#8d806780","#9dbc9b80")) +
  ylab("Biomass proportion") +
  xlab("multiplicative change in relative height") +
  theme_bw() +
  theme(legend.title = element_blank())