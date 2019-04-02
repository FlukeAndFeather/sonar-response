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

# Locomotion power (Pout) -----------------------------------------------------
fs_fun <- function(U, L) 1.5 * U / L
U_b_ms <- 1.5

CL_fun <- function(m) 1.46 + 0.0005 * m

Pout_fun <- function(u, l, m) {
  f_f <- fs_fun(u, l)
  f_b <- fs_fun(U_b_ms, l)
  (f_f - f_b) * CL_fun(m) * m
}

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
