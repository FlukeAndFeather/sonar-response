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

# prey_cdf_fun returns a CDF function paramaterized to prey data.
# Specifically, the probability function of energy per feeding event in
# kJ/gulp. Dave calculated the mean and standard deviation of the distribution
# of biomass density (kg/m3), log10 transformed. To get from biomass to energy:
# E       = b     * rho   * V
# kJ/gulp = kg/m3 * kJ/kg * m3/gulp
# given: mean(log10_b), sd(log_10b)
# log10(b) = log10(E / V / rho)
prey_cdf_fun <- function(mean_b, sd_b, V, rho) {
  function(E) {
    pnorm(log10(E / V / rho), 
          mean = mean_b,
          sd = sd_b)
  }
}
prey_quant_fun <- function(mean_b, sd_b, V, rho) {
  function(p, ...) {
    log10_b <- qnorm(p, 
                     mean = mean_b,
                     sd = sd_b)
    10 ^ log10_b * V * rho
  }
}

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

# Prey data
prey_data <- read_csv("data/Cetacea model output BOUT_EXTANT_v2.csv") %>% 
  select(-(X13:X25)) %>% 
  mutate(Species = recode(Species,
                          bonarensis = "bonaerensis",
                          Phocaena = "phocoena"),
         Genus = recode(Genus,
                        Physter = "Physeter",
                        Phocoaena = "Phocoena"),
         binomial = factor(paste(Genus, Species),
                           levels = binom_levels)) %>% 
  # data duplicated for multiple MR exponent values, choose one
  filter(`MR exponent` == 0.45)  

# Dave's rorqual data
# NOTE: engulfment volumes are from Potvin, not Kahane-Rapport
rorqual_prey_data <- read_csv("data/BaleenWhaleForagingDistBigKrill100Bins.csv") %>%
  select(1:10) %>% 
  slice(2:5) %>% 
  mutate(binomial = factor(recode(Species,
                                  bw = "Balaenoptera musculus",
                                  bp = "Balaenoptera physalus",
                                  mn = "Megaptera novaeangliae",
                                  bb = "Balaenoptera bonaerensis"),
                           levels = binom_levels),
         p_Ep = pmap(list(meanLog10biomass_log10kgm3,
                          stdLog10biomass_log10kgm3, 
                          meanVw_m3gulp,
                          KrillEnergy_kJkg),
                     prey_cdf_fun),
         q_Ep = pmap(list(meanLog10biomass_log10kgm3,
                          stdLog10biomass_log10kgm3, 
                          meanVw_m3gulp,
                          KrillEnergy_kJkg),
                     prey_quant_fun))


# For odontocetes, fit log-normal to Ep
prey_cdf_tbl <- prey_data %>% 
  group_by(binomial) %>% 
  group_map(function(data, key) {
    binom <- key[[1]]
    family <- data$Family[1]
    # ECDF data (Matt)
    emp_data <- uncount(data, Percent * 10)
    
    # For rorquals, use Dave's numbers
    if (family == "Balaenopteridae") {
      # Log-normal CDF data (Dave)
      q_Ep_fun <- filter(rorqual_prey_data, binomial == binom)$q_Ep[[1]]
      p_Ep_fun <- filter(rorqual_prey_data, binomial == binom)$p_Ep[[1]]
    } else {
      # For odontocetes, fit log-normal
      log10Ep_mean <- mean(log10(data$`Energy (kJ)`))
      log10Ep_sd <- sd(log10(data$`Energy (kJ)`))
      q_Ep_fun <- function(q) 10 ^ qnorm(q, 
                                         mean = log10Ep_mean, 
                                         sd = log10Ep_sd)
      p_Ep_fun <- function(Ep) pnorm(log10(Ep),
                                     mean = log10Ep_mean, 
                                     sd = log10Ep_sd)
    }
    
    lncdf_bounds <- q_Ep_fun(c(0.01, 0.99))
    lncdf_data <- tibble(Ep = seq(lncdf_bounds[1],
                                  lncdf_bounds[2],
                                  length.out = 100),
                         p_Ep = p_Ep_fun(Ep))
    
    plot <- ggplot() +
      # ECDF
      stat_ecdf(aes(`Energy (kJ)`),
                emp_data) +
      # Log-normal
      geom_line(aes(Ep, p_Ep),
                lncdf_data,
                linetype = "dashed") +
      labs(x = "Energy per feeding event (kJ)",
           y = "Cumulative probability",
           title = binom) +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Uncomment me to re-generate prey CDF plots
    # ggsave(sprintf("figs/prey_cdfs/%s.pdf", binom),
    #        plot,
    #        width = 9,
    #        height = 6)
    
    tibble(plot = list(plot),
           q_Ep_fun = list(q_Ep_fun),
           p_Ep_fun = list(p_Ep_fun))
  })

# Rorqual feeding rates
rorqual_data <- read_csv("data/lunge_rates_from_Paolo.csv") %>%
  filter(species != "be") %>% 
  mutate(duration_h = `deployment-time_secs` / 60 / 60,
         binomial = factor(recode(species,
                                  ba = "Balaenoptera bonaerensis",
                                  bp = "Balaenoptera physalus",
                                  bw = "Balaenoptera musculus",
                                  mn = "Megaptera novaeangliae"),
                           levels = binom_levels),
         lunge_h = total_lunges / duration_h)

# Odontocete feeding rates
odontocete_data <- read_csv("data/foragestats_combined_ko2.csv") %>%
  mutate(binomial = factor(str_replace(Species, "_", " "),
                           levels = binom_levels)) %>% 
  filter(taxa == "O")

# Feeding rates (rf) ----------------------------------------------------------
# Rorqual lunge rates
lunge_rate <- rorqual_data %>% 
  filter(duration_h >= 24) %>%
  group_by(binomial) %>% 
  summarize(N = n(),
            rf_h = round(mean(lunge_h), digits = 1)) %>% 
  mutate(Group = "Rorqual")

# Odontocete buzz rates
buzz_rate <- odontocete_data %>% 
  group_by(binomial) %>% 
  summarize(N = n(),
            rf_h = (total_buzz_count / total_duration_h) %>% 
              mean(na.rm = TRUE) %>% 
              round(digits = 1)) %>% 
  mutate(Group = "Odontocete")

# Feeding rates table
rf_tbl <- rbind(lunge_rate, buzz_rate)

# Once I get Danuta's buzz data I'll update this from a deployment-based
# distribution to hourly
rf_tbl2 <- rbind(
  transmute(rorqual_data,
            ID, 
            binomial, 
            rf_h = lunge_h),
  transmute(odontocete_data,
            ID, 
            binomial, 
            rf_h = total_buzz_count / total_duration_h)
) %>% 
  drop_na

# Feeding rate quantiles
buzz_tbl <- read_csv("data/deployment_info.csv") %>% 
  left_join(read_csv("data/all_buzzes.csv")) %>% 
  mutate(binomial = factor(str_replace(species, "_", " "),
                           levels = binom_levels)) %>% 
  group_by(deployment_ID) %>% 
  mutate(hour = floor(buzz_start / 3600)) %>% 
  group_by(binomial, deployment_ID, hour) %>% 
  summarize(rf_h = n())

buzz_q_tbl <- buzz_tbl %>% 
  group_by(binomial) %>% 
  group_map(function(data, key) {
    binom <- key[[1]]
    # ECDF plot
    plot <- ggplot(data, aes(rf_h)) + 
      stat_ecdf() +
      labs(x = "Hourly feeding rate",
           y = "Cumulative probability",
           title = binom) +
      theme_minimal()
    # ggsave(sprintf("figs/rf_cdfs/%s.pdf", binom),
    #        plot,
    #        width = 9,
    #        height = 6)
    # Quantile function
    q_buzz <- function(p, ...) {
      quantile(data$rf_h, probs = p)
    }
    
    tibble(plot = list(plot),
           q_rf_fun = list(q_buzz))
  })

# Prey energy (Ep) ------------------------------------------------------------
Ep_tbl <- prey_data %>%
  uncount(weights = Percent) %>% 
  group_by(binomial) %>% 
  summarize(Ep_kJ = mean(`Energy (kJ)`))

# From Dave, rorqual prey density is a log-normal distribution. Are odontocete
# prey densities purely empirical? Or can I estimate with log-normal as well?

# Consumption power (Pin) -----------------------------------------------------
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
# rf (Empirical)
# Quantiles of the empirical feed rate distribution
rf_q <- rbind(
  transmute(rorqual_data,
            ID, 
            binomial, 
            rf_h = lunge_h),
  transmute(odontocete_data,
            ID, 
            binomial, 
            rf_h = total_buzz_count / total_duration_h)
  ) %>% 
  drop_na %>% 
  group_by(binomial) %>% 
  group_map(~ tibble(q_rf = list(function(p, ...) quantile(.x$rf_h, p)))) %>% 
  ungroup

# Ep (log-normal, see prey_cdf_tbl)
# Ub (uniform, 0.5 - 2 m/s)
# ff_mult (uniform, 1/5x - 5x)
# CL (intercept and slope: uniform, 1/5x - 5x)

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
         ff_hz <- fs_fun(Uf_ms, L_m) * ff_mult
         fb_hz <- fs_fun(U_b_ms, L_m) * ff_mult
         Pout_W <- (ff_hz - fb_hz) * (CL_int + CL_slope * M_kg)
         
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
param <- c("rf_h", "Ep_kJ", "ff_mult", "CL_int", "CL_slope")
# List of parameter distribution functions
q <- list(rf_h = filter(rf_q, binomial == "Balaenoptera musculus")$q_rf[[1]],
          Ep_kJ = filter(prey_cdf_tbl, binomial == "Balaenoptera musculus")$q_Ep_fun[[1]],
          ff_mult = qunif,
          CL_int = qunif,
          CL_slope = qunif)
# List of distribution function parameters
q_arg <- list(rf_h = list(),
              Ep_kJ = list(),
              ff_mult = list(min = 0.2, max = 5),
              CL_int = list(min = 1.46 / 2, max = 1.46 * 2),
              CL_slope = list(min = 0.0005 / 2, max = 0.0005 * 2))
# Latin hypercube sample of parameter space
bw_LHS <- pse::LHS(Esonar_fun, param, 200, q, q_arg, nboot = 50)
# ECDF of model outputs
pse::plotecdf(bw_LHS)
# Scatter of model outputs w.r.t. parameters
pse::plotscatter(bw_LHS)

# Case studies ----------------------------------------------------------------
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
