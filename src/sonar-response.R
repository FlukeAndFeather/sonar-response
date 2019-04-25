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
Ein_fun <- function(td_min) {
  function(rf_h, Ep_kJ) {
    # Consumption power
    Pin_kJh <- rf_h * Ep_kJ
    Pin_W <- Pin_kJh / 3600
    
    # Consumption energy
    Ein_kJ <- Pin_kJh * td_min / 60
    
    Ein_kJ
  }
}
Eout_fun <- function(tf_min) {
  function(delta_ff, CL) {
    # Locomotion power
    Pout_W <- delta_ff * CL
    
    # Locomotion energy
    Eout_kJ <- Pout_W / 1000 * tf_min * 60
    
    Eout_kJ
  }
}
Esonar_fun <- function(td_min, tf_min) {
  function(data) {
    with(data, 
         {
           Ein_kJ <- Ein_fun(td_min)(rf_h, Ep_kJ)
           Eout_kJ <- Eout_fun(tf_min)(delta_ff, CL)
           
           Ein_kJ + Eout_kJ
         }
    )
  }
}

# Scenarios -------------------------------------------------------------------
## Behavioral responses
scenario_tbl <- tribble(
  ~binomial,                 ~t_d_min, ~t_f_min, ~U_f_ms,
  "Phocoena phocoena",       60,       60,       5,
  "Grampus griseus",         60,       60,       5,
  "Ziphius cavirostris",     60,       60,       5,
  "Physeter macrocephalus",  60,       60,       5,
  "Megaptera novaeangliae",  60,       60,       5,
  "Balaenoptera musculus",   60,       60,       5
) %>% 
  mutate(binomial = factor(binomial, levels = binom_levels)) %>% 
  ## Morphologies
  left_join(select(morphologies,
                   binomial,
                   Mass_kg,
                   Length_m),
            by = "binomial") %>%
  ## rf probabilities
  left_join(select(rf_tbl, binomial, q_rf_fun), by = "binomial") %>% 
  ## Ep probabilities
  left_join(select(Ep_tbl, binomial, q_Ep_fun), by = "binomial") 

Esonar_tbl <- scenario_tbl %>% 
  group_by(binomial) %>% 
  group_map(function(data, key) {
    param <- c("rf_h", "Ep_kJ", "delta_ff", "CL")
    q <- list(rf_h = data$q_rf_fun[[1]],
              Ep_kJ = data$q_Ep_fun[[1]],
              delta_ff = qgamma,
              CL = qgamma)
    
    # List of distribution function parameters
    ff_flight <- fs_fun(data$U_f_ms, data$Length_m)
    ff_basal <- fs_fun(U_b_ms, data$Length_m)
    delta_ff <- ff_flight - ff_basal
    ff_shape <- 4
    ff_scale <- delta_ff / ff_shape
    CL <- 1.46 + 0.0005 * data$Mass_kg
    CL_shape <- 4
    CL_scale <- CL / CL_shape
    q_arg <- list(rf_h = list(),
                  Ep_kJ = list(),
                  delta_ff = list(shape = ff_shape, scale = ff_scale),
                  CL = list(shape = CL_shape, scale = CL_scale))
    
    param_args <- tibble(delta_ff = delta_ff,
                         ff_shape = ff_shape,
                         ff_scale = ff_scale,
                         CL = CL,
                         CL_shape = CL_shape,
                         CL_scale = CL_scale)
    
    # Plots of delta_ff, CL distributions
    ggplot(data.frame(x = c(0, qgamma(0.99, 
                                      shape = ff_shape, 
                                      scale = ff_scale))), 
           aes(x)) +
      stat_function(fun = dgamma, 
                    args = q_arg[[3]]) +
      geom_vline(xintercept = delta_ff,
                 linetype = "dashed") +
      labs(x = "Change in fluking frequency (Hz)",
           y = "Probability density", 
           title = key$binomial) +
      theme_minimal()
    ggsave(sprintf("figs/ff_density/%s.pdf", key$binomial),
           width = 9,
           height = 6)
    ggplot(data.frame(x = c(0, qgamma(0.99, 
                                      shape = CL_shape, 
                                      scale = CL_scale))), 
           aes(x)) +
      stat_function(fun = dgamma, 
                    args = q_arg[[4]]) +
      geom_vline(xintercept = CL,
                 linetype = "dashed") +
      labs(x = "Locomotor cost (J/kg/stroke)",
           y = "density",
           title = key$binomial) +
      theme_minimal()
    ggsave(sprintf("figs/ff_density/%s.pdf", key$binomial),
           width = 9,
           height = 6)
    
    # Latin hypercube sample of parameter space
    model <- Esonar_fun(data$t_d_min,
                        data$t_f_min)
    esonar_LHS <- pse::LHS(model, param, 1e3, q, q_arg)
    
    sens_result <- esonar_LHS$data %>% 
      mutate(Esonar_kJ = esonar_LHS$res[,1,1],
             Eout_kJ = Eout_fun(data$t_f_min)(delta_ff, CL),
             Ein_kJ = Ein_fun(data$t_d_min)(rf_h, Ep_kJ),
             inout_ratio = Ein_kJ / Eout_kJ)
    
    # ECDF of model outputs
    esonar_ecdf <- ggplot(sens_result, aes(Esonar_kJ)) +
      stat_ecdf() +
      labs(x = "Energetic cost (kJ)",
           y = "Cumulative probability",
           title = key$binomial) +
      theme_minimal()
    
    # Scatter of model outputs w.r.t. parameters
    esonar_scatter <- sens_result %>% 
      gather(parameter, value, rf_h:CL) %>% 
      ggplot(aes(value, Esonar_kJ)) +
      geom_point(size = 0.5) +
      geom_smooth(method = "lm",
                  se = FALSE) +
      labs(x = "",
           y = "Energetic cost (kJ)",
           title = key$binomial) +
      facet_wrap(~ parameter, 
                 scales = "free_x",
                 strip.position = "bottom") +
      theme_minimal() +
      theme(strip.placement = "outside")
    
    # Save plots
    ggsave(sprintf("figs/esonar_ecdfs/%s.pdf", key$binomial),
           esonar_ecdf,
           width = 9,
           height = 6)
    ggsave(sprintf("figs/esonar_scatters/%s.pdf", key$binomial),
           esonar_scatter,
           width = 9,
           height = 6)
    
    # Linear model results
    esonar_linear <- sens_result %>%
      # Normalize values using z-scores
      mutate_at(vars(Esonar_kJ, rf_h, Ep_kJ, CL, delta_ff), 
                function(x) (x - mean(x)) / sd(x)) %>% 
      # Multiple regression
      lm(Esonar_kJ ~ rf_h + Ep_kJ + CL + delta_ff, data = .)
    
    # Extract coefficients, p-values, and confidence intervals
    esonar_coef <- coef(esonar_linear)
    esonar_pval <- summary(esonar_linear)$coefficients[,4]
    esonar_ci <- as_tibble(confint(esonar_linear, level = 0.95),
                           rownames = "param")
    colnames(esonar_ci)[2:3] <- c(ci_min, ci_max)
      
    # Combine and drop info for the intercept
    lm_results <- cbind(esonar_ci, esonar_coef, esonar_pval)[-1,]
    
    esonar_results <- summarize(sens_result,
                                mean_Esonar = mean(Esonar_kJ),
                                median_Esonar = median(Esonar_kJ),
                                iqr_Esonar = IQR(Esonar_kJ),
                                median_inout = median(inout_ratio),
                                inout_25 = quantile(inout_ratio, 0.25),
                                inout_75 = quantile(inout_ratio, 0.75))
    
    cbind(lm_results, esonar_results, param_args)
  })

# A plot with the normalized linear model coefficients
coef_data <- Esonar_tbl %>% 
  mutate(param = factor(param, 
                        levels = c("CL",
                                   "delta_ff",
                                   "Ep_kJ",
                                   "rf_h"),
                        labels = c("C[L]",
                                   "Delta*f[f]",
                                   "E[p]",
                                   "r[f]")))
ggplot(coef_data, aes(x = binomial, y = esonar_coef, group = binomial)) +
  geom_hline(aes(yintercept = mean_coef),
             summarize(group_by(coef_data, param), 
                       mean_coef = mean(esonar_coef)),
             linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                width = 0.2, 
                position = position_dodge(width = 0.4)) +
  facet_grid(~ param,
             labeller = label_parsed,
             switch = "x") +
  labs(x = "Parameter",
       y = "Normalized coefficient") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside")
ggsave("figs/extreme_flight_coefs.pdf",
       width = 9,
       height = 6)

# Ratio of energy in to energy out
ggplot(Esonar_tbl, aes(binomial, median_inout)) +
  geom_point() +
  geom_errorbar(aes(ymin = inout_25+1, ymax = inout_75),
                width = 0.2) +
  scale_y_log10(labels = log_labels) +
  scale_x_discrete(labels = function(lbl) str_replace(lbl, " ", "\n")) +
  labs(y = expression(E["in"]:E[out])) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Case studies ----------------------------------------------------------------
casestudy_tbl <- tribble(
  ~binomial,                 ~t_d_min, ~t_f_min, ~U_f_ms,
  "Ziphius cavirostris",     6 * 60,   30,       4.5,
  "Balaenoptera musculus",   1 * 60,   5,        3
) %>% 
  mutate(binomial = factor(binomial, levels = binom_levels)) %>% 
  ## Morphologies
  left_join(select(morphologies,
                   binomial,
                   Mass_kg,
                   Length_m),
            by = "binomial") %>%
  ## rf probabilities
  left_join(select(rf_tbl, binomial, q_rf_fun), by = "binomial") %>% 
  ## Ep probabilities
  left_join(select(Ep_tbl, binomial, q_Ep_fun), by = "binomial") 

Esonar_tbl2 <- casestudy_tbl %>% 
  group_by(binomial) %>% 
  group_map(function(data, key) {
    param <- c("rf_h", "Ep_kJ", "delta_ff", "CL")
    q <- list(rf_h = data$q_rf_fun[[1]],
              Ep_kJ = data$q_Ep_fun[[1]],
              delta_ff = qgamma,
              CL = qgamma)
    
    # List of distribution function parameters
    ff_flight <- fs_fun(data$U_f_ms, data$Length_m)
    ff_basal <- fs_fun(U_b_ms, data$Length_m)
    delta_ff <- ff_flight - ff_basal
    ff_shape <- 4
    ff_scale <- delta_ff / ff_shape
    CL <- 1.46 + 0.0005 * data$Mass_kg
    CL_shape <- 4
    CL_scale <- CL / CL_shape
    q_arg <- list(rf_h = list(),
                  Ep_kJ = list(),
                  delta_ff = list(shape = ff_shape, scale = ff_scale),
                  CL = list(shape = CL_shape, scale = CL_scale))
    
    param_args <- tibble(delta_ff = delta_ff,
                         ff_shape = ff_shape,
                         ff_scale = ff_scale,
                         CL = CL,
                         CL_shape = CL_shape,
                         CL_scale = CL_scale)
    
    # Plots of delta_ff, CL distributions
    ggplot(data.frame(x = c(0, qgamma(0.99, 
                                      shape = ff_shape, 
                                      scale = ff_scale))), 
           aes(x)) +
      stat_function(fun = dgamma, 
                    args = q_arg[[3]]) +
      geom_vline(xintercept = delta_ff,
                 linetype = "dashed") +
      labs(x = "Change in fluking frequency (Hz)",
           y = "Probability density", 
           title = key$binomial) +
      theme_minimal()
    ggsave(sprintf("figs/ff_density/%s.pdf", key$binomial),
           width = 9,
           height = 6)
    ggplot(data.frame(x = c(0, qgamma(0.99, 
                                      shape = CL_shape, 
                                      scale = CL_scale))), 
           aes(x)) +
      stat_function(fun = dgamma, 
                    args = q_arg[[4]]) +
      geom_vline(xintercept = CL,
                 linetype = "dashed") +
      labs(x = "Locomotor cost (J/kg/stroke)",
           y = "density",
           title = key$binomial) +
      theme_minimal()
    ggsave(sprintf("figs/ff_density/%s.pdf", key$binomial),
           width = 9,
           height = 6)
    
    # Latin hypercube sample of parameter space
    model <- Esonar_fun(data$t_d_min,
                        data$t_f_min)
    esonar_LHS <- pse::LHS(model, param, 1e3, q, q_arg)
    
    sens_result <- esonar_LHS$data %>% 
      mutate(Esonar_kJ = esonar_LHS$res[,1,1],
             Eout_kJ = Eout_fun(data$t_f_min)(delta_ff, CL),
             Ein_kJ = Ein_fun(data$t_d_min)(rf_h, Ep_kJ),
             inout_ratio = Ein_kJ / Eout_kJ)
    
    # ECDF of model outputs
    esonar_ecdf <- ggplot(sens_result, aes(Esonar_kJ)) +
      stat_ecdf() +
      labs(x = "Energetic cost (kJ)",
           y = "Cumulative probability",
           title = key$binomial) +
      theme_minimal()
    
    # Scatter of model outputs w.r.t. parameters
    esonar_scatter <- sens_result %>% 
      gather(parameter, value, rf_h:CL) %>% 
      ggplot(aes(value, Esonar_kJ)) +
      geom_point(size = 0.5) +
      geom_smooth(method = "lm",
                  se = FALSE) +
      labs(x = "",
           y = "Energetic cost (kJ)",
           title = key$binomial) +
      facet_wrap(~ parameter, 
                 scales = "free_x",
                 strip.position = "bottom") +
      theme_minimal() +
      theme(strip.placement = "outside")
    
    # Save plots
    ggsave(sprintf("figs/esonar_ecdfs/%s.pdf", key$binomial),
           esonar_ecdf,
           width = 9,
           height = 6)
    ggsave(sprintf("figs/esonar_scatters/%s.pdf", key$binomial),
           esonar_scatter,
           width = 9,
           height = 6)
    
    # Linear model results
    esonar_linear <- sens_result %>%
      # Normalize values using z-scores
      mutate_at(vars(Esonar_kJ, rf_h, Ep_kJ, CL, delta_ff), 
                function(x) (x - mean(x)) / sd(x)) %>% 
      # Multiple regression
      lm(Esonar_kJ ~ rf_h + Ep_kJ + CL + delta_ff, data = .)
    
    # Extract coefficients, p-values, and confidence intervals
    esonar_coef <- coef(esonar_linear)
    esonar_pval <- summary(esonar_linear)$coefficients[,4]
    esonar_ci <- as_tibble(confint(esonar_linear, level = 0.95),
                           rownames = "param")
    colnames(esonar_ci)[2:3] <- c(ci_min, ci_max)
    
    # Combine and drop info for the intercept
    lm_results <- cbind(esonar_ci, esonar_coef, esonar_pval)[-1,]
    
    esonar_results <- summarize(sens_result,
                                mean_Esonar = mean(Esonar_kJ),
                                median_Esonar = median(Esonar_kJ),
                                iqr_Esonar = IQR(Esonar_kJ),
                                median_inout = median(inout_ratio),
                                inout_25 = quantile(inout_ratio, 0.25),
                                inout_75 = quantile(inout_ratio, 0.75))
    
    cbind(lm_results, esonar_results, param_args)
  })


# A plot with the normalized linear model coefficients
coef_data2 <- Esonar_tbl2 %>% 
  mutate(param = factor(param, 
                        levels = c("CL",
                                   "delta_ff",
                                   "Ep_kJ",
                                   "rf_h"),
                        labels = c("C[L]",
                                   "Delta*f[f]",
                                   "E[p]",
                                   "r[f]")))
ggplot(coef_data, aes(x = binomial, y = esonar_coef, group = binomial)) +
  geom_hline(aes(yintercept = mean_coef),
             summarize(group_by(coef_data, param), 
                       mean_coef = mean(esonar_coef)),
             linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                width = 0.2, 
                position = position_dodge(width = 0.4)) +
  facet_grid(~ param,
             labeller = label_parsed,
             switch = "x") +
  labs(x = "Parameter",
       y = "Normalized coefficient") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside")

ggsave("figs/extreme_flight_coefs.pdf",
       width = 9,
       height = 6)