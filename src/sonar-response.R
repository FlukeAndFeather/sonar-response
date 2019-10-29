library(ggrepel)
library(ggsci)
library(lmodel2)
library(scales)
library(tidyverse)

# Utility functions ----------------------------------------------------------

# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom = function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

# Labels for logarithmic scales
log_labels <- trans_format("log10", math_format(10^.x))

cbf_palette <- c("Phocoenidae and Delphinidae" = rgb(0, 114, 178, maxColorValue = 255), # blue (descent)
                 "Physeteridae and Ziphiidae" = rgb(213, 94, 0, maxColorValue = 255), # vermillion (ascent)
                 "Balaenopteridae" = rgb(0, 158, 115, maxColorValue = 255)) # bluish green (surface)

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
  # binomial is a factor ordered by species length 
  arrange(Length_m) %>% 
  mutate(binomial = factor(Species, 
                           levels = unique(Species)),
         abbr = str_split(binomial, " ") %>% 
           map_chr(~ paste0(str_sub(.x[1], 1, 1), str_sub(.x[2], 1, 1))),
         abbr = case_when(binomial == "Globicephala melas" ~ "Gme",
                          binomial == "Globicephala macrorhynchus" ~ "Gma",
                          TRUE ~ abbr)) %>% 
  filter(!binomial %in% c("Orcinus orca", "Berardius bairdii"))
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
  filter(!binomial %in% c("Orcinus orca", "Berardius bairdii")) %>% 
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
# Ep/m figure
prey_tbl %>% 
  filter(!binomial %in% c("Orcinus orca", "Berardius bairdii")) %>% 
  mutate(binomial = factor(binomial, levels = binom_levels),
         grouping = case_when(Family %in% c("Phocoenidae", "Delphinidae") ~ "Phocoenidae and Delphinidae",
                              Family %in% c("Physeteridae", "Ziphiidae") ~ "Physeteridae and Ziphiidae",
                              Family == "Balaenopteridae" ~ "Balaenopteridae")) %>% 
  left_join(select(morphologies, binomial, abbr, Mass_kg),
            by = "binomial") %>% 
  mutate(meanEp_kJkg = meanEp_kJ / Mass_kg,
         firstqEp_kJkg = firstqEp_kJ / Mass_kg,
         thirdqEp_kJkg = thirdqEp_kJ / Mass_kg) %>% 
  ggplot(aes(Mass_kg, meanEp_kJkg, color = grouping, shape = grouping)) +
  geom_pointrange(aes(ymin = firstqEp_kJkg, ymax = thirdqEp_kJkg),
                  fatten = 3,
                  size = 0.75) +
  geom_text_repel(aes(label = abbr),
                  size = 3) +
  scale_x_log10(labels = log_labels) +
  scale_y_log10() +
  scale_color_manual(values = cbf_palette) +
  labs(x = "Body mass (kg)",
       y = expression("Mass-specific " * E[p] ~ (kJ ~ kg ^ -1))) +
  theme_classic(base_size = 12) +
  theme(axis.title = element_text(size = 10),
        legend.position = "none")
ggsave("figs/Ep2.pdf",
       width = 80,
       height = 65,
       units = "mm",
       dpi = 600)

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
  filter(binomial != "Orcinus orca") %>% 
  ggplot(aes(binomial, mean_rf, color = grouping, shape = grouping)) +
  geom_pointrange(aes(ymin = firstq_rf, ymax = thirdq_rf),
                  fatten = 3,
                  size = 0.75) +
  scale_x_discrete(labels = morphologies$abbr) +
  scale_color_manual(values = cbf_palette) +
  labs(y = expression("Feeding rate " * (hr^-1))) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 45,
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    margin = margin(0, 0, 0, 0)),
        legend.position = "none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5))
ggsave("figs/rf.pdf",
       width = 80,
       height = 65,
       units = "mm",
       dpi = 600)

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
                     labels = function(x) x / 10^5) +
  labs(x = expression(italic(E[p]) ~ (10^5 ~ kJ))) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8))
bw_rf <- filter(rf_tbl, binomial == "Balaenoptera musculus")
rf_inset <- ggplot(tibble(x = bw_rf$q_rf_fun[[1]](seq(0, 1, length.out = 1000))),
                   aes(x)) + 
  geom_histogram(binwidth = 1, 
                 boundary = 0, 
                 fill = "light gray", 
                 color = "black",
                 size = 0.2) +
  labs(x = expression(italic(r[f]) ~ ("hr"^{-1}))) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8))
bw_Pc <- filter(Pin_tbl, binomial == "Balaenoptera musculus") %>% 
  group_by(binomial) %>% 
  group_modify(function(data, key) {
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
                     name = expression(italic(frac(dE[a], dt)) ~ ("kJ" ~ "hr"^{-1}))) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8))
Pc_plot +
  annotation_custom(ggplotGrob(ep_inset), 
                    xmin = 0.7e7, xmax = 2.45e7,
                    ymin = 90, ymax = 360) +
  annotation_custom(ggplotGrob(rf_inset), 
                    xmin = 2.5e7, xmax = 4.25e7,
                    ymin = 90, ymax = 360)
# Dimensions and DPI for Conservation Letters
# https://authorservices.wiley.com/asset/photos/electronic_artwork_guidelines.pdf
ggsave("figs/Pc.pdf", 
       width = 80, 
       height = 65,
       units = "mm",
       dpi = 600)

# Pc table
Pc_tbl <- rf_tbl %>% 
  left_join(Ep_tbl, by = "binomial") %>% 
  left_join(select(morphologies, binomial, Mass_kg)) %>% 
  group_by(binomial) %>% 
  group_modify(function(data, key) {
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
  group_modify(function(data, key) {
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
                              Family == "Balaenopteridae" ~ "Balaenopteridae") %>% 
           factor(levels = c("Phocoenidae and Delphinidae",
                             "Physeteridae and Ziphiidae",
                             "Balaenopteridae")))

# For presentations -----------------------------------------------------------
consumption_coef <- filter(coef_data, scenario == "consumption")
# Consumption, all species
ggplot(consumption_coef, aes(x = param, y = esonar_coef)) +
  geom_boxplot() +
  coord_flip() +
  scale_x_discrete(labels = parse(text = levels(consumption_coef$param))) +
  scale_y_continuous(limits = c(-0.25, 1.0),
                     breaks = seq(-0.25, 1.0, by = 0.25)) +
  labs(y = "Sensitivity") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())
ggsave("figs/consumption_allsp.pdf",
       width = 4.5,
       height = 3)
# Consumption, by grouping
consumption_grouped <- ggplot(
  consumption_coef, 
  aes(x = param, y = esonar_coef, color = grouping)
) +
  geom_boxplot() +
  coord_flip() +
  scale_x_discrete(labels = parse(text = levels(consumption_coef$param))) +
  scale_y_continuous(limits = c(-0.25, 1.0),
                     breaks = seq(-0.25, 1.0, by = 0.25)) +
  scale_color_manual(values = cbf_palette) +
  labs(y = "Sensitivity") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
consumption_grouped
ggsave("figs/consumption_groups.pdf",
       width = 4.5,
       height = 3)

flight_coef <- filter(coef_data, scenario == "flight")
# Flight, all species
ggplot(flight_coef, aes(x = param, y = esonar_coef)) +
  geom_boxplot() +
  coord_flip() +
  scale_x_discrete(labels = parse(text = levels(flight_coef$param))) +
  scale_y_continuous(limits = c(-0.25, 1.0),
                     breaks = seq(-0.25, 1.0, by = 0.25)) +
  labs(y = "Sensitivity") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())
ggsave("figs/flight_allsp.pdf",
       width = 4.5,
       height = 3)
# Flight, by grouping
flight_grouped <- ggplot(flight_coef, 
                        aes(x = param, y = esonar_coef, color = grouping)) +
  geom_boxplot() +
  coord_flip() +
  scale_x_discrete(labels = parse(text = levels(flight_coef$param))) +
  scale_y_continuous(limits = c(-0.25, 1.0),
                     breaks = seq(-0.25, 1.0, by = 0.25)) +
  scale_color_manual(values = cbf_palette) +
  labs(y = "Sensitivity") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
flight_grouped
ggsave("figs/flight_groups.pdf",
       width = 4.5,
       height = 3)

# Figure for paper
consumption_grouped + 
  theme(axis.ticks = element_line(),
        axis.text.x = element_text()) +
  flight_grouped + 
  theme(axis.ticks = element_line(),
        axis.text.x = element_text()) +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A")
ggsave("figs/sensitivity.pdf",
       width = 180,
       height = 80,
       units = "mm",
       dpi = 600)
# Ratio of in/out
coef_data %>% 
  mutate(inout = ifelse(param %in% c("r[f]", "E[p]"), "in", "out")) %>% 
  group_by(inout, scenario) %>% 
  summarize(mean_coef = mean(esonar_coef))

ggplot(consumption_coef, aes(x = binomial, y = esonar_coef)) +
  # geom_vline(aes(xintercept = mean_coef),
  #            summarize(group_by(consumption_coef, param), 
  #                      mean_coef = mean(esonar_coef)),
  #            linetype = "dashed") +
  geom_boxplot() +
  facet_grid(param ~ .,
             labeller = label_parsed,
             switch = "y") +
  scale_x_continuous(limits = c(-0.25, 1.0),
                     breaks = seq(-0.25, 1.0, by = 0.25)) +
  labs(x = "Sensitivity") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")

# A table of the ratio of energy in to energy out
Esonar_tbl %>% 
  group_by(binomial, scenario) %>% 
  slice(1) %>% 
  ungroup %>% 
  select(binomial, scenario, mean_inout, se_inout) %>% 
  write_csv("figs/flight_inout.csv")

# Energy in vs energy out
inout_fun <- function(t_d_min, t_f_min, U_f_ms, binom) {
  morph_row <- filter(morphologies, binomial == binom)
  M_kg <- morph_row$Mass_kg
  L_m <- morph_row$Length_m
  # Parameter distributions
  rf_fun <- filter(rf_tbl, binomial == binom)$q_rf_fun[[1]]
  Ep_fun <- function(p) {
    row <- filter(prey_tbl, binomial == binom)
    qlnorm(p, meanlog = row$meanlnEp_lnkJ, sdlog = row$sdlnEp_lnkJ)
  }
  ff_fun <- function(p) {
    if (U_f_ms == 0) {
      0
    } else {
      ff_flight <- fs_fun(U_f_ms, L_m)
      ff_basal <- fs_fun(U_b_ms, L_m)
      delta_ff <- ff_flight - ff_basal
      ff_shape <- 4
      ff_scale <- delta_ff / ff_shape
      qgamma(p, ff_shape, scale = ff_scale)
    }
  }
  CL_fun <- function(p) {
    CL <- 1.46 + 0.0005 * M_kg
    CL_shape <- 4
    CL_scale <- CL / CL_shape
    qgamma(p, CL_shape, scale = CL_scale)
  }
  
  rf_q <- runif(1)
  rf_h <- rf_fun(rf_q)
  Ep_q <- runif(1)
  Ep_kJ <- Ep_fun(Ep_q)
  Ein <- Ein_fun(t_d_min)(rf_h, Ep_kJ)
  ff_q <- runif(1)
  delta_ff <- ff_fun(ff_q)
  CL_q <- runif(1)
  CL <- CL_fun(CL_q)
  Eout <- Eout_fun(t_f_min, M_kg)(delta_ff, CL)
  Esonar <- Ein + Eout
  tibble(rf_q, rf_h, Ep_q, Ep_kJ, ff_q, delta_ff, CL_q, CL, Ein, Eout, Esonar)
}

Einout <- tribble(
  ~scenario,        ~t_d_min, ~t_f_min, ~U_f_ms,
  "cessation",      60,       0,        0,
  "mild_flight",    0,        5,        2.5,
  "strong_flight",  0,        15,       3.5,
  "extreme_flight", 0,        30,       5
) %>% 
  mutate(scenario = factor(scenario,
                           levels = c("cessation", 
                                      "mild_flight",
                                      "strong_flight",
                                      "extreme_flight"),
                           labels = c("Cessation only",
                                      "Mild flight",
                                      "Strong flight",
                                      "Extreme flight"))) %>% 
  # All but Oo and Bb
  crossing(binomial = binom_levels[c(-6, -9)]) %>% 
  mutate(binomial = factor(binomial, levels = binom_levels)) %>% 
  group_by(binomial, scenario, t_d_min, t_f_min, U_f_ms, binomial) %>% 
  group_modify(function(data, key) {
    map_dfr(1:500, ~ inout_fun(key$t_d_min, key$t_f_min, key$U_f_ms, key$binomial))
  })

Einout %>% 
  ungroup %>% 
  left_join(morphologies, by = "binomial") %>% 
  mutate(grouping = case_when(Family %in% c("Phocoenidae", "Delphinidae") ~ "Phocoenidae and Delphinidae",
                              Family %in% c("Physeteridae", "Ziphiidae") ~ "Physeteridae and Ziphiidae",
                              Family == "Balaenopteridae" ~ "Balaenopteridae"),
         scenario = fct_relabel(scenario, ~str_replace(.x, " ", "\n"))) %>% 
  ggplot(aes(scenario, Esonar, color = grouping)) +
  geom_boxplot() +
  scale_color_brewer(palette = "Dark2") +
  labs(y = "Energy cost (kJ)") +
  facet_wrap(~ binomial,
             scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave("figs/inout.pdf",
       width = 9,
       height = 6)

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
  # All but Oo and Bb
  crossing(binomial = binom_levels[c(-6, -9)],
           beta = c(2.5, 3, 5),
           bmr = factor(c("Kleiber", "Maresh")))

thr_calculator <- function(U_f_ms, tf_min, binom, beta, bmr) {
  # Morphologies
  morph_row <- filter(morphologies, binomial == binom)
  if (nrow(morph_row) != 1) stop("Imprecise binomial")
  L_m = morph_row$Length_m
  M_kg = morph_row$Mass_kg
  
  
  # Parameter distributions
  rf_fun <- filter(rf_tbl, binomial == binom)$q_rf_fun[[1]]
  Ep_fun <- function(p) {
    row <- filter(prey_tbl, binomial == binom)
    qlnorm(p, meanlog = row$meanlnEp_lnkJ, sdlog = row$sdlnEp_lnkJ)
  }
  ff_fun <- function(p) {
    if (U_f_ms == 0) {
      0
    } else {
      ff_flight <- fs_fun(U_f_ms, L_m)
      ff_basal <- fs_fun(U_b_ms, L_m)
      delta_ff <- ff_flight - ff_basal
      ff_shape <- 4
      ff_scale <- delta_ff / ff_shape
      qgamma(p, ff_shape, scale = ff_scale)
    }
  }
  CL_fun <- function(p) {
    CL <- 1.46 + 0.0005 * M_kg
    CL_shape <- 4
    CL_scale <- CL / CL_shape
    qgamma(p, CL_shape, scale = CL_scale)
  }
  
  # Daily FMR 
  daily_bmr = if (bmr == "Kleiber") {
    293.1 * M_kg ^ 0.75
  } else if (bmr == "Maresh") {
    581 * M_kg ^ 0.68
  } else {
    stop("Unspecified BMR calculation")
  }
  daily_FMR <- beta * daily_bmr
  
  thr_fun <- function() {
    Esonar <- 0
    model <- Esonar_fun(60, tf_min, M_kg)
    td_min <- 0
    #print(sprintf("%s: tf_min=%.1f, beta=%.1f, bmr=%s", binom, tf_min, beta, bmr))
    while (Esonar < daily_FMR) {
      rf_h <- rf_fun(runif(1))
      Ep_kJ <- Ep_fun(runif(1))
      delta_ff <- ff_fun(runif(1))
      CL <- CL_fun(runif(1))
      
      result <- model(tibble(rf_h, Ep_kJ, delta_ff, CL))
      if (Esonar + result < daily_FMR) {
        td_min <- td_min + 60
      } else {
        td_min <- td_min + 60 * (daily_FMR - Esonar) / result
      }
      Esonar <- Esonar + result
    }
    td_min
  }
  
  td_min <- map_dbl(1:1e3, ~ thr_fun())
  
  tibble(mean_thr = mean(td_min),
         med_thr = median(td_min),
         firstq_thr = quantile(td_min, 0.25),
         thirdq_thr = quantile(td_min, 0.75))
}

thr_tbl <- scenario_tbl %>% 
  mutate(binomial = factor(binomial, levels = binom_levels)) %>% 
  group_by_all %>% 
  group_modify(function(data, key) {
    with(key, thr_calculator(U_f_ms, t_f_min, binomial, beta, bmr))
  }) %>% 
  arrange(binomial, t_f_min)

# filter(thr_tbl, beta == 3, bmr == "Maresh") %>% 
filter(thr_tbl, beta == 3) %>% 
  ggplot(aes(binomial, 
             mean_thr, 
             color = scenario,
             shape = bmr)) +
  geom_pointrange(aes(ymin = firstq_thr,
                      ymax = thirdq_thr),
                  fatten = 2,
                  position = position_dodge(width = 0.6)) +
  scale_x_discrete(labels = morphologies$abbr) +
  scale_y_continuous(breaks = c(30, 60, 60 * 4, 60 * 12, 60*24, 60*48, 60*72),
                     minor_breaks = NULL,
                     labels = c("30 min", "1 hour", "4 hours", "12 hours", "1 day", "2 days", "3 days"),
                     trans = "log2") +
  scale_color_brewer(palette = "RdYlBu", direction = -1) +
  labs(y = "Feeding cessation") +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, -10))
ggsave("figs/critical_threshold.pdf",
       width = 180,
       height = 120,
       units = "mm",
       dpi = 600)
 
 # Threshold table
 thr_tbl %>% 
   filter(beta == 3) %>% 
   write_csv("data/output/thresholds.csv")
   
   
