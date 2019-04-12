library(tidyverse)
library(scales)
library(lmodel2)
library(ggsci)

# Utility functions ----------------------------------------------------------

# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom = function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

# Labels for logarithmic scales
log_labels <- trans_format("log10", math_format(10^.x))

# Data ------------------------------------------------------------------------
# Morphological data
morphologies <- read_csv("data/foragestats_combined_ko2.csv") %>% 
  mutate(Species = str_replace(Species, "_", " ")) %>% 
  group_by(Species) %>% 
  summarize(Length_m = first(Body_length_m),
            Mass_kg = first(Body_mass_kg)) %>% 
  mutate(Clade = ifelse(str_detect(Species, ".*ptera.*"),
                        "Mysticete",
                        "Odontocete"),
         Family = recode(Species,
                         `Balaenoptera bonaerensis` = "Balaenopteridae",
                         `Balaenoptera musculus` = "Balaenopteridae",
                         `Balaenoptera physalus` = "Balaenopteridae",
                         `Berardius bairdii` = "Ziphiidae",
                         `Globicephala macrorhynchus` = "Delphinidae",
                         `Globicephala melas` = "Delphinidae",
                         `Grampus griseus` = "Delphinidae",
                         `Megaptera novaeangliae` = "Balaenopteridae",
                         `Mesoplodon densirostris` = "Ziphiidae",
                         `Orcinus orca` = "Delphinidae",
                         `Phocoena phocoena` = "Phocoenidae",
                         `Physeter macrocephalus` = "Physeteridae",
                         `Ziphius cavirostris` = "Ziphiidae")) %>% 
  # binomial is a factor ordered by species length with an abbreviated label
  arrange(Length_m) %>% 
  mutate(binomial = factor(Species, 
                           levels = unique(Species)))
binom_levels <- levels(morphologies$binomial)

# Prey
load("data/prey_cdf_tbl.RData")

# Feeding rates
load("data/buzz_rf.RData")
load("data/lunge_rf.RData")

# Consumption power (Pin) -----------------------------------------------------
Ep_tbl <- prey_cdf_tbl %>% 
  select(binomial, Ep_kJ = mean_Ep, q_Ep_fun) %>% 
  ungroup %>% 
  mutate(binomial = factor(binomial, levels = binom_levels))
rf_tbl <- rbind(buzz_rf, lunge_rf) %>% 
  select(binomial, rf_h = mean_rf, q_rf_fun) %>% 
  ungroup %>% 
  mutate(binomial = factor(binomial, levels = binom_levels))
Pin_tbl <- inner_join(Ep_tbl, rf_tbl, by = "binomial") %>% 
  mutate(Pin_kJ_h = Ep_kJ * rf_h) %>% 
  left_join(select(morphologies, binomial, Family), by = "binomial")

# Locomotion power (Pout) -----------------------------------------------------
fs_fun <- function(U, L, La = 0.2, St = 0.3) {
  # St = A * f / U
  # A = L * La
  # St = L * La * f / U
  # f = St * U / L / La
  St * U / L / La
}
U_b_ms <- 1.5

CL_fun <- function(m) 1.46 + 0.0005 * m

Pout_fun <- function(u, l, m) {
  f_f <- fs_fun(u, l)
  f_b <- fs_fun(U_b_ms, l)
  (f_f - f_b) * CL_fun(m) * m
}

# Sensitivity -----------------------------------------------------------------
# Parameter CDFs
# rf (Empirical, see lunge_times.R and buzz_times.R)
# Ep (log-normal, see prey_cdf.R)
# Ub (uniform, 0.5 - 2 m/s)
# Both fluke frequency and locomotor cost multipliers are used as 10^mult,
# shifting the base value by an order of magnitude in both directions
# ff_mult (uniform, -10 - 10)
# CL_mult (uniform, -10 - 10)

# Vectorized function for calculating Esonar for sensitivity analysis
Esonar_fun <- function(data) {
  M_kg <- 93000
  L_m <- 25.20
  td_min <- 60
  tf_min <- 5
  Uf_ms <- 2.5
  with(data, 
       {
         # Consumption and locomotion power
         Pin_kJh <- rf_h * Ep_kJ
         Pin_W <- Pin_kJh / 3600
         ff_hz <- fs_fun(Uf_ms, L_m) * 10 ^ ff_mult
         fb_hz <- fs_fun(U_b_ms, L_m)
         Pout_W <- (ff_hz - fb_hz) * 10 ^ CL_mult * (1.46 + 0.0005 * M_kg)
         
         # Consumption and locomotion energy (based on behavioral responses)
         Ein_kJ <- Pin_kJh * td_min / 60
         Eout_kJ <- Pout_W / 1000 * tf_min * 60
         
         # Esonar in kJ
         Ein_kJ + Eout_kJ
       }
  )
}

# Sensitivity analysis using pse package
# List of model parameters
param <- c("rf_h", "Ep_kJ", "ff_mult", "CL_mult")
# List of parameter distribution functions
bw_rf_q <- filter(rf_tbl, binomial == "Balaenoptera musculus")$q_rf_fun[[1]]
bw_Ep_q <- filter(Ep_tbl, binomial == "Balaenoptera musculus")$q_Ep_fun[[1]]
q <- list(rf_h = bw_rf_q,
          Ep_kJ = bw_Ep_q,
          ff_mult = qunif,
          CL_mult = qunif)
# List of distribution function parameters
q_arg <- list(rf_h = list(),
              Ep_kJ = list(),
              ff_mult = list(min = -1, max = 1),
              CL_mult = list(min = -1, max = 1))
# Latin hypercube sample of parameter space
bw_LHS <- pse::LHS(Esonar_fun, param, 500, q, q_arg, nboot = 50)
sens_result <- cbind(bw_LHS$data, bw_LHS$res)
colnames(sens_result)[length(colnames(sens_result))] <- "Esonar"
# ECDF of model outputs
ggplot(sens_result, aes(Esonar)) +
  stat_ecdf() +
  labs(x = "Energetic cost (kJ)",
       y = "Cumulative probability") +
  theme_minimal()
# Scatter of model outputs w.r.t. parameters
pse::plotscatter(bw_LHS)
sens_result %>% 
  gather(parameter, value, rf_h:CL_mult) %>% 
  ggplot(aes(value, Esonar)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm",
              se = FALSE) +
  labs(x = "",
       y = "Energetic cost (kJ)") +
  facet_wrap(~ parameter, 
             scales = "free_x",
             strip.position = "bottom") +
  theme_minimal() +
  theme(strip.placement = "outside")

# Scenarios -------------------------------------------------------------------
cases_tbl <- tibble(binomial = factor(c("Mesoplodon densirostris",
                                        "Ziphius cavirostris",
                                        "Balaenoptera bonaerensis",
                                        "Balaenoptera musculus"),
                                      levels = binom_levels),
                    t_d_min = c(6*60,
                                6*60,
                                2.5*60,
                                60),
                    t_f_min = c(30,
                                30,
                                60,
                                5),
                    U_f_ms = c(4.5,
                               4.5,
                               3.5,
                               2.5),
                    Reference = c("DeRuiter et al. 2013 Fig 1.",
                                  "DeRuiter et al. 2013 Fig 1.",
                                  "Kvadsheim et al. 2017 Fig 2.",
                                  "Southall et al. 2019 bw11\\_219b"))

esonar_tbl <- cases_tbl %>% 
  left_join(select(morphologies,
                   binomial,
                   Mass_kg,
                   Length_m),
            by = "binomial") %>%
  # P_in
  left_join(Pin_tbl, by = "binomial") %>%
  # P_out
  mutate(f_f = fs_fun(U_f_ms, Length_m),
         f_b = fs_fun(U_b_ms, Length_m),
         C_L = CL_fun(Mass_kg),
         Pout_W = Pout_fun(U_f_ms, Length_m, Mass_kg)) %>% 
  # E_sonar
  # The units are a little wonky. Intake is in units of kJ/hour and output
  # is in J/s (W).
  mutate(E_in_kJ = Pin_kJ_h * t_d_min / 60,
         E_out_kJ = Pout_W / 1000 * t_f_min * 60,
         E_sonar_kJ = E_in_kJ + E_out_kJ,
         BMR_kJ_day = 293.1 * Mass_kg ^ 0.75,
         E_BMR = E_sonar_kJ / BMR_kJ_day,
         E_m = E_sonar_kJ / Mass_kg)
