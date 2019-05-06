library(tidyverse)

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
         # Change base from 10 to e
         meanlnbiomass_lnkgm3 = log(10 ^ meanLog10biomass_log10kgm3),
         sdlnbiomass_lnkgm3 = log(10 ^ (stdLog10biomass_log10kgm3 ^ 2)) ^ 0.5,
         # Multiply biomass by volume and energy density to get energy
         meanlnenergy_lnkJ = meanlnbiomass_lnkgm3 + log(meanVw_m3gulp * KrillEnergy_kJkg),
         sdlnenergy_lnkJ = sdlnbiomass_lnkgm3)

# For odontocetes, fit log-normal to Ep
prey_tbl <- prey_data %>% 
  group_by(Family, binomial) %>% 
  group_map(function(data, key) {
    binom <- key$binomial[1]
    family <- key$Family[1]
    
    # For rorquals, use Dave's numbers
    if (family == "Balaenopteridae") {
      # Log-normal params (Dave)
      species_prey_data <- filter(rorqual_prey_data, binomial == binom)
      meanlnEp_lnkJ <- species_prey_data$meanlnenergy_lnkJ
      sdlnEp_lnkJ <- species_prey_data$sdlnenergy_lnkJ
    } else {
      # For odontocetes, fit log-normal
      lnenergy_lnkj <- log(data$`Energy (kJ)`)
      energy_perc <- data$Percent
      meanlnEp_lnkJ <- Hmisc::wtd.mean(lnenergy_lnkj, energy_perc)
      sdlnEp_lnkJ <- sqrt(Hmisc::wtd.var(lnenergy_lnkj, energy_perc))
    }
    
    # Function for plotting Ep distribution
    Ep_plot <- function(statfun, ylab) {
      ggplot(tibble(Ep = qlnorm(c(0, 0.99), meanlnEp_lnkJ, sdlnEp_lnkJ)), 
             aes(Ep)) +
        stat_function(fun = statfun, args = list(meanlnEp_lnkJ, sdlnEp_lnkJ)) +
        geom_vline(xintercept = exp(meanlnEp_lnkJ), linetype = "dashed") +
        labs(x = "Energy per feeding event (kJ)",
             y = ylab,
             title = binom) +
        theme_minimal()
    } 
    
    # Plot density function
    Ep_plot(dlnorm, "Probability density")
    ggsave(sprintf("figs/prey_density/%s.pdf", binom),
           width = 9,
           height = 6)
    
    # Plot cumulative distribution function
    Ep_plot(plnorm, "Cumulative probability")
    ggsave(sprintf("figs/prey_cdfs/%s.pdf", binom),
           width = 9,
           height = 6)
    
    tibble(meanEp_kJ = exp(meanlnEp_lnkJ),
           firstqEp_kJ = qlnorm(c(0.25), meanlnEp_lnkJ, sdlnEp_lnkJ),
           thirdqEp_kJ = qlnorm(c(0.75), meanlnEp_lnkJ, sdlnEp_lnkJ),
           iqrEp_kJ = thirdqEp_kJ - firstqEp_kJ,
           meanlnEp_lnkJ = meanlnEp_lnkJ,
           sdlnEp_lnkJ = sdlnEp_lnkJ)
  }) %>% 
  ungroup

save(prey_tbl, file = "data/prey_tbl.RData")
