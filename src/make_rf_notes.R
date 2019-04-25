library(tidyverse)

good_lunges <- c("good", "good_dives", "good dives", "dives good")

# Only good quality lunges, not in Chile, no sonar exposure
rf_notes <- read_csv("data/feeding_rate_notes-7-30-18.csv") %>% 
  drop_na %>% 
  filter(site != "Chile",
         lunge_quality %in% good_lunges,
         sonar_exp == "no")

write_csv(rf_notes, "data/rf_notes.csv")
