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
load("data/prey_tbl.RData")

# Feeding rates
load("data/buzz_rf.RData")
load("data/Md_buzz_rf.RData")
load("data/lunge_rf.RData")

# Consumption power (Pin) -----------------------------------------------------
Ep_tbl <- prey_tbl %>% 
  select(binomial, 
         meanEp_kJ, 
         meanlnEp_lnkJ,
         sdlnEp_lnkJ,
         firstqEp_kJ,
         thirdqEp_kJ) %>% 
  ungroup %>% 
  mutate(binomial = factor(binomial, levels = binom_levels))

# Ep figure
prey_tbl %>% 
  mutate(binomial = factor(binomial, levels = binom_levels),
         grouping = case_when(Family %in% c("Phocoenidae", "Delphinidae") ~ "Phocoenidae and Delphinidae",
                              Family %in% c("Physeteridae", "Ziphiidae") ~ "Physeteridae and Ziphiidae",
                              Family == "Balaenopteridae" ~ "Balaenopteridae")) %>% 
  ggplot(aes(binomial, meanEp_kJ, color = grouping)) +
  geom_point() +
  geom_errorbar(aes(ymin = firstqEp_kJ, ymax = thirdqEp_kJ),
                width = 0.4) +
  scale_x_discrete(labels = function(lbl) str_replace(lbl, " ", "\n")) +
  scale_y_log10(labels = log_labels) +
  scale_color_aaas() +
  labs(y = "Energy per feeding event (kJ)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave("figs/Ep.pdf",
       width = 9,
       height = 6)

rf_tbl <- bind_rows(buzz_rf, lunge_rf, Md_buzz_rf) %>% 
  select(binomial, 
         rf_h = mean_rf, 
         firstq_rf,
         thirdq_rf,
         q_rf_fun) %>% 
  ungroup %>% 
  mutate(binomial = factor(binomial, levels = binom_levels))

# rf figure
bind_rows(buzz_rf, lunge_rf, Md_buzz_rf) %>% 
  select(binomial, mean_rf, firstq_rf, thirdq_rf) %>% 
  left_join(select(prey_tbl, binomial, Family), by = "binomial") %>% 
  mutate(binomial = factor(binomial, levels = binom_levels),
         grouping = case_when(Family %in% c("Phocoenidae", "Delphinidae") ~ "Phocoenidae and Delphinidae",
                              Family %in% c("Physeteridae", "Ziphiidae") ~ "Physeteridae and Ziphiidae",
                              Family == "Balaenopteridae" ~ "Balaenopteridae")) %>% 
  ggplot(aes(binomial, mean_rf, color = grouping)) +
  geom_point() +
  geom_errorbar(aes(ymin = firstq_rf, ymax = thirdq_rf),
                width = 0.4) +
  scale_x_discrete(labels = function(lbl) str_replace(lbl, " ", "\n")) +
  scale_color_aaas() +
  labs(y = "Feeding rate (events / hour)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave("figs/rf.pdf",
       width = 9,
       height = 6)

Pin_tbl <- inner_join(Ep_tbl, rf_tbl, by = "binomial") %>% 
  mutate(Pin_kJ_h = meanEp_kJ * rf_h) %>% 
  left_join(select(morphologies, binomial, Family, Mass_kg), by = "binomial")

# Mass-specific consumption power
sample_Pin <- function(rf_q, meanlnEp, sdlnEp, n = 1e3) {
  pse::LHS(function(data) data$rf * data$Ep,
           factors = c("rf", "Ep"),
           N = n,
           q = list(rf = rf_q,
                    Ep = qlnorm),
           q.arg = list(rf = list(),
                        Ep = list(meanlog = meanlnEp,
                                  sdlog = sdlnEp)),
           res.names = "Pin_kJhr")$res[,,1]
}

# Figure 1, Pc is bimodally distributed
# Using Bw as example
bw_ep <- filter(prey_tbl, binomial == "Balaenoptera musculus")
ep_inset <- ggplot(tibble(x = qlnorm(c(0.001, 0.99),
                                     meanlog = bw_ep$meanlnEp_lnkJ,
                                     sdlog = bw_ep$sdlnEp_lnkJ)),
                   aes(x)) +
  stat_function(fun = dlnorm, 
                args = list(meanlog = bw_ep$meanlnEp_lnkJ,
                            sdlog = bw_ep$sdlnEp_lnkJ)) +
  scale_x_continuous(breaks = seq(0, 1e6, by = 250e3),
                     labels = c(0, 
                                expression(2.5 %*% 10^5), 
                                expression(5 %*% 10^5), 
                                expression(7.5 %*% 10^5), 
                                expression(1 %*% 10^6))) +
  labs(x = expression(italic(E[p]) ~~ (kJ))) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
bw_rf <- filter(rf_tbl, binomial == "Balaenoptera musculus")
rf_inset <- ggplot(tibble(x = bw_rf$q_rf_fun[[1]](seq(0, 1, length.out = 1000))),
                   aes(x)) + 
  geom_histogram(binwidth = 1, 
                 boundary = 0, 
                 fill = "light gray", 
                 color = "black") +
  labs(x = expression(italic(r[f]) ~~ ("hr"^{-1}))) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
bw_Pc <- filter(Pin_tbl, binomial == "Balaenoptera musculus") %>% 
  group_by(binomial) %>% 
  group_map(function(data, key) {
    with(data,
         tibble(Pc = sample_Pin(q_rf_fun[[1]],
                                meanlnEp_lnkJ[1],
                                sdlnEp_lnkJ[1])))
  }) %>% 
  ungroup
Pc_plot <- ggplot(bw_Pc, aes(Pc)) +
  geom_histogram(bins = 30,
                 boundary = 0, 
                 fill = "light gray", 
                 color = "black") +
  geom_vline(aes(xintercept = mean(Pc)),
             linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 4e7, by = 1e7),
                     limits = c(0, 4e7),
                     labels = c(0,
                                expression(1 %*% 10^7),
                                expression(2 %*% 10^7),
                                expression(3 %*% 10^7),
                                expression(4 %*% 10^7)),
                     name = expression(italic(P[c]) ~~ ("kJ" %.% "hr"^{-1}))) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(1,2,1,1)/5, "in"))
Pc_plot +
  annotation_custom(ggplotGrob(ep_inset), 
                    xmin = 1e7, xmax = 2.25e7,
                    ymin = 160, ymax = 360) +
  annotation_custom(ggplotGrob(rf_inset), 
                    xmin = 2.35e7, xmax = 3.6e7,
                    ymin = 160, ymax = 360)
ggsave("figs/Pc.pdf", width = 9, height = 6)

# Pc table
Pc_tbl <- rf_tbl %>% 
  left_join(Ep_tbl, by = "binomial") %>% 
  left_join(select(morphologies, binomial, Mass_kg)) %>% 
  group_by(binomial) %>% 
  group_map(function(data, key) {
    Pc <- sample_Pin(data$q_rf_fun[[1]],
                     data$meanlnEp_lnkJ,
                     data$sdlnEp_lnkJ)
    tibble(meanEp_kJ = data$meanEp_kJ,
           meanEp_kJkg = meanEp_kJ / data$Mass_kg,
           iqr1Ep_kJ = data$firstqEp_kJ,
           iqr3Ep_kJ = data$thirdqEp_kJ,
           meanrf_h = data$rf_h,
           meanrf_hkg = meanrf_h / data$Mass_kg,
           iqr1rf_h = data$firstq_rf,
           iqr3rf_h = data$thirdq_rf,
           meanPc_kJh = mean(Pc),
           meanPc_kJhkg = meanPc_kJh / data$Mass_kg,
           iqr1Pc_kJh = quantile(Pc, 0.25),
           iqr3Pc_kJh = quantile(Pc, 0.75))
  })
write_csv(Pc_tbl, "data/output/Pc.csv")

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
Eout_fun <- function(tf_min, m) {
  function(delta_ff, CL) {
    # Locomotion power
    Pout_W <- delta_ff * CL * m
    
    # Locomotion energy
    Eout_kJ <- Pout_W / 1000 * tf_min * 60
    
    Eout_kJ
  }
}
Esonar_fun <- function(td_min, tf_min, m) {
  function(data) {
    with(data, 
         {
           Ein_kJ <- Ein_fun(td_min)(rf_h, Ep_kJ)
           Eout_kJ <- Eout_fun(tf_min, m)(delta_ff, CL)
           
           Ein_kJ + Eout_kJ
         }
    )
  }
}

# Scenarios -------------------------------------------------------------------
# Extreme flight
# Behavioral responses
scenario_tbl <- 
  tribble(~scenario,     ~t_d_min, ~t_f_min, ~U_f_ms,
          "flight",      60,       30,       5,
          "consumption", 240,      10,       3.5) %>% 
  crossing(select(morphologies,
                  binomial,
                  Mass_kg,
                  Length_m) %>% 
             filter(!binomial %in% c("Orcinus orca",
                                     "Berardius bairdii"))) %>% 
  # rf probabilities
  left_join(select(rf_tbl, binomial, q_rf_fun), by = "binomial") %>% 
  # Ep probabilities
  left_join(select(Ep_tbl, binomial, meanlnEp_lnkJ, sdlnEp_lnkJ), 
            by = "binomial") 

Esonar_tbl <- scenario_tbl %>% 
  group_by(scenario, binomial) %>% 
  group_map(function(data, key) {
    param <- c("rf_h", "Ep_kJ", "delta_ff", "CL")
    q <- list(rf_h = data$q_rf_fun[[1]],
              Ep_kJ = qlnorm,
              delta_ff = qgamma,
              CL = qgamma)
    
    # List of distribution function parameters
    ff_flight <- fs_fun(data$U_f_ms, data$Length_m)
    ff_basal <- fs_fun(U_b_ms, data$Length_m)
    delta_ff <- ff_flight - ff_basal
    ff_shape <- 4
    ff_scale <- delta_ff / ff_shape
    # squared to get units of J/stroke
    CL <- 1.46 + 0.0005 * data$Mass_kg
    CL_shape <- 4
    CL_scale <- CL / CL_shape
    q_arg <- list(rf_h = list(),
                  Ep_kJ = list(meanlog = data$meanlnEp_lnkJ, 
                               sdlog = data$sdlnEp_lnkJ),
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
           title = key$binomial,
           caption = sprintf("U_b = %.1f m/s, U_f = %.1f m/s, f_b = %.2f Hz, f_f = %.2f Hz",
                             U_b_ms,
                             data$U_f_ms,
                             ff_basal,
                             ff_flight)) +
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
      labs(x = "Locomotor cost (J/stroke)",
           y = "density",
           title = key$binomial) +
      theme_minimal()
    ggsave(sprintf("figs/CL_density/%s.pdf", key$binomial),
           width = 9,
           height = 6)
    
    # Latin hypercube sample of parameter space
    model <- Esonar_fun(data$t_d_min,
                        data$t_f_min,
                        data$Mass_kg)
    tryCatch(esonar_LHS <- pse::LHS(model, param, 1e2, q, q_arg),
             error = function(e) browser())
    
    sens_result <- esonar_LHS$data %>% 
      mutate(Esonar_kJ = esonar_LHS$res[,1,1],
             Eout_kJ = Eout_fun(data$t_f_min, data$Mass_kg)(delta_ff, CL),
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
    ggsave(sprintf("figs/esonar_ecdfs/%s_%s.pdf", key$scenario, key$binomial),
           esonar_ecdf,
           width = 9,
           height = 6)
    ggsave(sprintf("figs/esonar_scatters/%s_%s.pdf", key$scenario, key$binomial),
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
    colnames(esonar_ci)[2:3] <- c("ci_min", "ci_max")
      
    # Combine and drop info for the intercept
    lm_results <- cbind(esonar_ci, esonar_coef, esonar_pval)[-1,]
    
    esonar_results <- summarize(sens_result,
                                mean_Esonar = mean(Esonar_kJ),
                                median_Esonar = median(Esonar_kJ),
                                iqr_Esonar = IQR(Esonar_kJ),
                                median_inout = median(inout_ratio),
                                inout_25 = quantile(inout_ratio, 0.25),
                                inout_75 = quantile(inout_ratio, 0.75),
                                mean_inout = mean(inout_ratio),
                                se_inout = sd(inout_ratio)/sqrt(n()))
    
     cbind(lm_results, esonar_results)
  })

# A plot with the normalized linear model coefficients
coef_data <- Esonar_tbl %>% 
  ungroup %>% 
  mutate(param = factor(param, 
                        levels = c("CL",
                                   "delta_ff",
                                   "Ep_kJ",
                                   "rf_h"),
                        labels = c("C[L]",
                                   "Delta*f[f]",
                                   "E[p]",
                                   "r[f]"))) %>% 
  left_join(select(morphologies, binomial, Family),
            by = "binomial") %>% 
  mutate(grouping = case_when(Family %in% c("Phocoenidae", "Delphinidae") ~ "Phocoenidae and Delphinidae",
                              Family %in% c("Physeteridae", "Ziphiidae") ~ "Physeteridae and Ziphiidae",
                              Family == "Balaenopteridae" ~ "Balaenopteridae"))

coef_plots <- coef_data %>% 
  group_by(scenario) %>% 
  group_split() %>% 
  map(~ ggplot(.x, aes(x = binomial, y = esonar_coef, color = grouping)) +
        geom_hline(aes(yintercept = mean_coef),
                   summarize(group_by(.x, param), 
                             mean_coef = mean(esonar_coef)),
                   linetype = "dashed") +
        geom_point(position = position_dodge(width = 0.4)) +
        geom_errorbar(aes(ymin = ci_min, ymax = ci_max), 
                      width = 0.2, 
                      position = position_dodge(width = 0.4)) +
        scale_color_brewer(palette = "Dark2") +
        facet_grid(~ param,
                   labeller = label_parsed,
                   switch = "x") +
        scale_y_continuous(limits = c(-0.25, 1.0),
                           breaks = seq(-0.25, 1.0, by = 0.25)) +
        labs(y = "Normalized coefficient") +
        theme_classic() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              legend.position = "none",
              legend.title = element_blank(),
              strip.background = element_blank(),
              strip.placement = "outside"))
coef_plots[[1]] + coef_plots[[2]] + patchwork::plot_layout(nrow = 1)
ggsave("figs/flight_coefs.pdf",
       width = 9,
       height = 6)
      
# A table of the ratio of energy in to energy out
Esonar_tbl %>% 
  group_by(binomial, scenario) %>% 
  slice(1) %>% 
  ungroup %>% 
  select(binomial, scenario, mean_inout, se_inout) %>% 
  write_csv("figs/flight_inout.csv")

# Critical thresholds ---------------------------------------------------------
scenario_tbl <- tribble(
  ~scenario,        ~t_f_min, ~U_f_ms,
  "no_flight",      0,        0,
  "mild_flight",    5,        2.5,
  "strong_flight",  15,       3.5,
  "extreme_flight", 30,       5
  ) %>% 
  mutate(scenario = factor(scenario,
                           levels = c("no_flight", 
                                      "mild_flight",
                                      "strong_flight",
                                      "extreme_flight"),
                           labels = c("No flight",
                                      "Mild flight",
                                      "Strong flight",
                                      "Extreme flight"))) %>% 
  crossing(binomial = c("Phocoena phocoena", 
                        "Grampus griseus", 
                        "Ziphius cavirostris", 
                        "Physeter macrocephalus", 
                        "Megaptera novaeangliae", 
                        "Balaenoptera musculus"),
           beta = c(2.5, 3, 5),
           bmr = factor(c("Kleiber", "Maresh")))

thr_factory <- function(U_f_ms, tf_min, binom, beta, bmr) {
  function(td_min) {
    # Morphologies
    morph_row <- filter(morphologies, binomial == binom)
    if (nrow(morph_row) != 1) stop("Imprecise binomial")
    L_m = morph_row$Length_m
    M_kg = morph_row$Mass_kg
    
    # Eout parameters 
    ff_flight <- fs_fun(U_f_ms, L_m)
    ff_basal <- fs_fun(U_b_ms, L_m)
    delta_ff <- ff_flight - ff_basal
    CL <- 1.46 + 0.0005 * M_kg
    
    # Ein parameters
    rf_row <- filter(rf_tbl, binomial == binom)
    if (nrow(rf_row) != 1) stop("Imprecise binomial")
    rf_h <- rf_row$q_rf_fun[[1]](0.5)
    
    Ep_row <- filter(Ep_tbl, binomial == binom)
    if (nrow(Ep_row) != 1) stop("Imprecise binomial")
    Ep_kJ <- Ep_row$q_Ep_fun[[1]](0.5)
    
    # Model result
    # Eout
    Pout_W <- delta_ff * CL * M_kg
    Eout_kJ <- Pout_W / 1000 * tf_min * 60
    # Ein
    Pin_kJh <- rf_h * Ep_kJ
    Pin_W <- Pin_kJh / 3600
    Ein_kJ <- Pin_kJh * td_min / 60
    
    Esonar <- Eout_kJ + Ein_kJ
    
    # Daily FMR 
    daily_bmr = if (bmr == "Kleiber") {
      293.1 * M_kg ^ 0.75
    } else if (bmr == "Maresh") {
      581 * M_kg ^ 0.68
    } else {
      stop("Unspecified BMR calculation")
    }
    daily_FMR <- beta * daily_bmr
    
    abs(Esonar - daily_FMR)
  }
}

find_thr <- function(thr_fun) {
  optimize(thr_fun, lower = 0, upper = 31 * 24 * 60)$minimum
}

thr_tbl <- scenario_tbl %>% 
  mutate(binomial = factor(binomial, levels = binom_levels),
         thr_fun = pmap(list(U_f_ms, t_f_min, binomial, beta, bmr), 
                        thr_factory),
         td_min_thr = map_dbl(thr_fun, find_thr),
         td_hr_thr = td_min_thr / 60,
         td_day_thr = td_hr_thr / 24) %>% 
  arrange(binomial, t_f_min)

filter(thr_tbl, beta == 3) %>% 
  ggplot(aes(binomial, 
             td_min_thr, 
             color = scenario,
             shape = bmr)) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_discrete(labels = function(lbl) str_replace(lbl, " ", "\n")) +
  scale_y_continuous(breaks = c(60, 60*24, 60*24*7, 60*24*7*2),
                     labels = c("Hour", "Day", "Week", "2 Weeks"),
                     trans = "log2") +
  # remove green to avoid color-blind unfriendly palette
  scale_color_manual(values = pal[-3]) +
  labs(y = "Time displaced from feeding") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("figs/critical_threshold.pdf",
       width = 9,
       height = 6)
