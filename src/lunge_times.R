library(zeallot)
library(tidyverse)
library(R.matlab)

dn_to_posix <- function(dn) {
  as.POSIXct((dn - 719529) * 86400, origin = "1970-01-01", tz = "UTC")
}

find_data <- function(id) {
  # Look for CATS
  data_dir <- dir("/Volumes/COPYCATS/CATS/tag_data", 
                  pattern = id, 
                  full.names = TRUE)
  list(prh_path = dir(data_dir,
                      pattern = "10Hzprh",
                      full.names = TRUE),
       lunge_path = dir(file.path(data_dir, "lunges"),
                        pattern = "lunges.*mat",
                        full.names = TRUE))
}

lunge_times <- function(data, key) {
  id <- key[[1]]
  data_paths <- find_data(id)
  c(prh, lunges) %<-% map(data_paths, readMat)
  
  tibble(lunge_dt = dn_to_posix(prh$DN[lunges$LungeI]))
}

deployments <- read_csv("data/lunge_rates_from_Paolo.csv") %>% 
  filter(lunge_quality %in% c("good", "good dives", "good_dives"), 
         sonar_exp == "none",
         species %in% c("ba", "bp", "bw", "mn"))

deployments %>% 
  group_by(ID) %>% 
  group_map(lunge_times)
