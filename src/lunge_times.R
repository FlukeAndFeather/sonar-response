library(zeallot)
library(tidyverse)
library(R.matlab)

# Convert MATLAB DN's to POSIXct
dn_to_posix <- function(dn) {
  as.POSIXct((dn - 719529) * 86400, origin = "1970-01-01", tz = "UTC")
}

# Look up PRH & lunge file paths on the CATS drive
# This will look different on Windows
find_data <- function(id) {
  # Look for CATS
  data_dir <- dir("/Volumes/COPYCATSdat/CATS/tag_data", 
                  pattern = id, 
                  full.names = TRUE)
  list(prh_path = dir(data_dir,
                      pattern = "10Hzprh",
                      full.names = TRUE),
       lunge_path = dir(file.path(data_dir, "lunges"),
                        pattern = "lunges.*mat",
                        full.names = TRUE))
}

# Find times of lunges. To be used with group_map. Takes an ID, finds the
# associated PRH and lunge files, converts times to POSIX, and returns a tibble
# with lunge times.
lunge_times <- function(data, key) {
  id <- key$ID
  data_paths <- find_data(id)
  tryCatch({
    c(prh, lunges) %<-% map(data_paths, readMat)
    tibble(lunge_dt = dn_to_posix(prh$DN[lunges$LungeI]))
  },
  error = function(e) tibble(lunge_dt = NA))
}

# Finds tag on and tag off times. To be used with group_map. Takes an ID, loads
# the PRH, finds the tag on/off times, and converts to POSIX. Returns a tibble
# with two columns (tagon and tagoff).
tag_times <- function(data, key) {
  id <- key$ID
  data_paths <- find_data(id)
  tryCatch({
    prh <- readMat(data_paths$prh_path)
    tagonoff <- dn_to_posix(range(prh$DN[prh$tagon == 1]))
    tibble(tagon = tagonoff[1],
           tagoff = tagonoff[2])
  }, error = function(e) {
    tibble(tagon = NA,
           tagoff = NA)
  })
}

# Deployment metadata. Starts with Paolo's lunge rate file. Filters down to
# good quality dives and the four rorqual species of interest. Removes sonar
# exposure deployments. Gets tag on/off times.
deployments <- read_csv("data/lunge_rates_from_Paolo.csv") %>% 
  filter(lunge_quality %in% c("good", "good dives", "good_dives"), 
         sonar_exp == "none",
         species %in% c("ba", "bp", "bw", "mn")) %>% 
  group_by(species, ID) %>% 
  group_map(tag_times)

# Looks up the lunge times for all the deployments.
lunge_tbl <- deployments %>% 
  group_by_all %>% 
  group_map(lunge_times)
 
# Final product. A tibble with columns for whale ID, species, tag on/off time, 
# hour of deployment, and feeding rate for that hour. Fills in missing hours.
# Drops the leading part of the deployment less than an hour. For example, 
# drop the first 0.5 hours of a 2.5 hour deployment. 
# Interesting note: approximately 25% of hours have feeding rate of 0 across 
# species.
hourly_lunges <- lunge_tbl %>% 
  drop_na() %>% 
  mutate(begin = tagon + (as.numeric(tagoff - tagon, units = "secs") %% 3600),
         hour = floor(as.numeric(lunge_dt - begin, units = "hours"))) %>% 
  filter(hour >= 0) %>% 
  group_by_at(vars(-lunge_dt)) %>% 
  summarize(rf_h = n()) %>%
  group_by_at(vars(-hour, -rf_h)) %>% 
  group_map(function(data, key) {
    dep_hrs <- floor(as.numeric(key$tagoff - key$tagon, units = "hours"))
    tibble(hour = seq_len(dep_hrs)) %>% 
      left_join(data, by = "hour") %>% 
      replace_na(list(rf_h = 0))
  }) %>% 
  ungroup %>% 
  mutate(binomial = recode(species, 
                           ba = "Balaenoptera acutorostrata",
                           bp = "Balaenoptera physalus",
                           bw = "Balaenoptera musculus",
                           mn = "Megaptera novaeangliae"))

# Save ECDF plots
hourly_lunges %>% 
  group_by(binomial) %>% 
  group_walk(function(data, key) {
    binom = key$binomial
    plot <- ggplot(data) +
      stat_ecdf(aes(rf_h)) +
      labs(x = "Hourly feeding rate",
           y = "Cumulative probability",
           title = binom) +
      theme_minimal()
    ggsave(sprintf("figs/rf_cdfs/%s.pdf", binom),
           plot,
           width = 9,
           height = 6)
  })

lunge_rf <- hourly_lunges %>% 
  # change ba to bb
  mutate(binomial = if_else(binomial == "Balaenoptera acutorostrata",
                            "Balaenoptera bonaerensis",
                            binomial)) %>% 
  group_by(binomial) %>% 
  group_map(function(data, key) {
    q_lunge <- function(p, ...) {
      quantile(data$rf_h, probs = p)
    }
    
    tibble(mean_rf = mean(data$rf_h),
           sd_rf = sd(data$rf_h),
           median_rf = median(data$rf_h),
           firstq_rf = q_lunge(0.25),
           thirdq_rf = q_lunge(0.75),
           q_rf_fun = list(q_lunge))
  }) %>% 
  ungroup

# Save raw and summarized data
save(hourly_lunges, file = "data/hourly_lunges.RData")
save(lunge_rf, file = "data/lunge_rf.RData")
