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

# Prey energy (Ep) ------------------------------------------------------------
Ep_tbl <- prey_data %>%
  uncount(weights = Percent) %>% 
  group_by(binomial) %>% 
  summarize(Ep_kJ = mean(`Energy (kJ)`))

# Consumption power (Pin) -----------------------------------------------------
Pin_tbl <- inner_join(Ep_tbl, rf_tbl, by = "binomial") %>% 
  mutate(Pin_kJ_h = Ep_kJ * rf_h) %>% 
  left_join(select(morphologies, binomial, Family), by = "binomial")

