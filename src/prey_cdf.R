library(tidyverse)

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

# Prey data
prey_data <- read_csv("data/Cetacea model output BOUT_EXTANT_v2.csv") %>% 
  select(-(X13:X25)) %>% 
  mutate(Species = recode(Species,
                          bonarensis = "bonaerensis",
                          Phocaena = "phocoena"),
         Genus = recode(Genus,
                        Physter = "Physeter",
                        Phocoaena = "Phocoena"),
         binomial = paste(Genus, Species)) %>% 
  # data duplicated for multiple MR exponent values, choose one
  filter(`MR exponent` == 0.45)  

# Dave's rorqual data
# NOTE: engulfment volumes are from Potvin, not Kahane-Rapport
rorqual_prey_data <- read_csv("data/BaleenWhaleForagingDistBigKrill100Bins.csv") %>%
  select(1:10) %>% 
  slice(2:5) %>% 
  mutate(binomial = recode(Species,
                           bw = "Balaenoptera musculus",
                           bp = "Balaenoptera physalus",
                           mn = "Megaptera novaeangliae",
                           bb = "Balaenoptera bonaerensis"),
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
      energy_kj <- log10(data$`Energy (kJ)`)
      energy_perc <- data$Percent
      log10Ep_mean <- Hmisc::wtd.mean(energy_kj, energy_perc)
      log10Ep_sd <- Hmisc::wtd.var(energy_kj, energy_perc) ^ 0.5
      q_Ep_fun <- function(p) 10 ^ qnorm(p, 
                                         mean = log10Ep_mean, 
                                         sd = log10Ep_sd)
      p_Ep_fun <- function(Ep) pnorm(log10(Ep),
                                     mean = log10Ep_mean, 
                                     sd = log10Ep_sd)
    }
    
    # Integrate Ep*p(Ep)*dEp to get mean(Ep)
    mean_Ep_tbl <- tibble(Ep = seq(q_Ep_fun(0.01), 
                                   q_Ep_fun(0.99), 
                                   length.out = 1e4),
                          p_Ep = p_Ep_fun(Ep),
                          d_Ep = lead(p_Ep) - p_Ep) %>% 
      summarize(sum(Ep * d_Ep, na.rm = TRUE))
    mean_Ep <- mean_Ep_tbl[[1]]

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
      # Vertical line at the mean
      geom_vline(xintercept = mean_Ep,
                 linetype = "dashed") +
      labs(x = "Energy per feeding event (kJ)",
           y = "Cumulative probability",
           title = binom) +
      theme_minimal() +
      theme(legend.position = "none")
    
    ggsave(sprintf("figs/prey_cdfs/%s.pdf", binom),
           plot,
           width = 9,
           height = 6)
    
    tibble(plot = list(plot),
           q_Ep_fun = list(q_Ep_fun),
           p_Ep_fun = list(p_Ep_fun),
           mean_Ep = mean_Ep,
           med_Ep = q_Ep_fun(0.5))
  })

save(prey_cdf_tbl, file = "data/prey_cdf_tbl.RData")
