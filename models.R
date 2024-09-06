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
library(tidybayes)
library(modelr)
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
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 12),
                      backend = "cmdstanr",
                      file = "fits/m.Bet_alle.rsl",
                      seed = 123)
summary(m.Bet_alle.rsl)
pp_check(m.Bet_alle.rsl, ndraws = 100)
conditional_effects(m.Bet_alle.rsl,
                    categorical = TRUE)



within.sp %>%
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
  theme(legend.title = element_blank())



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




m.all.rsl_hhm2 <- brm(bf(RSL ~ (rel.ht  + leaf_habit + myc_group)^2 + (1 + rel.ht|SppName/Study)), 
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





library(ape)

phylo = read.tree("phyliptree1.phy")

A = ape::vcv.phylo(phylo)
dimnames(A) = list(unique(within.sp$SppName),
                   unique(within.sp$SppName))

within.sp$phylo = within.sp$SppName

m.all.rsl_hhm_ph <- brm(bf(RSL ~ (rel.ht  + leaf_habit + myc_group)^2 + 
                             (1|gr(phylo, cov = A))
                           + (1|SppName)
), 
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
file = "fits/m.all.rsl_hhm_ph",
seed = 123)
summary(m.all.rsl_hhm2)