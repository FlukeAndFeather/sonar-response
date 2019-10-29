library(zeallot)
library(tidyverse)
library(R.matlab)
library(readxl)
library(lubridate)

# Convert MATLAB DN's to POSIXct
dn_to_posix <- function(dn) {
  as.POSIXct((dn - 719529) * 86400, origin = "1970-01-01", tz = "UTC")
}

# Look up PRH & lunge file paths on the CATS drive
# This will look different on Windows
find_cats_data <- function(id) {
  # Look for CATS
  cats_dir <- dir("/Volumes/COPYCATSdat/CATS/tag_data", 
                  pattern = id, 
                  full.names = TRUE)
  list(prh_path = dir(cats_dir,
                      pattern = "10Hzprh",
                      full.names = TRUE),
       lunge_path = dir(file.path(cats_dir, "lunges"),
                        pattern = "lunges.*mat",
                        full.names = TRUE))
}
find_dtag_data <- function(id) {
  # Look for dtag
  dtag_dir <- "/Volumes/COPYCATSdat/DTAG"
  list(prh_path = dir(file.path(dtag_dir, "prh"),
                      pattern = sprintf("%s.*", id),
                      full.names = TRUE),
       lunge_path = dir(file.path(dtag_dir, "lunges"),
                        pattern = sprintf("%s.*", id),
                        full.names = TRUE))
}

# Find times of lunges. To be used with group_map. Takes an ID, finds the
# associated PRH and lunge files, converts times to POSIX, and returns a tibble
# with lunge times.d
lunge_times_cats <- function(data, key) {
  id <- key$ID
  data_paths <- find_cats_data(id)
  tryCatch({
    c(prh, lunges) %<-% map(data_paths, readMat)
    tibble(lunge_dt = dn_to_posix(prh$DN[lunges$LungeI]))
  },
  error = function(e) tibble(lunge_dt = NA))
}
lunge_times_dtag <- function(data, key) {
  id <- key$ID
  tagon <- key$`Tag on`
  data_paths <- find_dtag_data(id)
  tryCatch({
    c(prh, lunges) %<-% map(data_paths, readMat)
    fs <- prh$fs[1]
    tibble(lunge_dt = (lunges$LungeI - 1) * fs + tagon)
  },
  error = function(e) tibble(lunge_dt = NA))
}

# Finds tag on and tag off times. To be used with group_map. Takes an ID, loads
# the PRH, finds the tag on/off times, and converts to POSIX. Returns a tibble
# with two columns (tagon and tagoff).
tag_times_cats <- function(data, key) {
  id <- key$ID
  data_paths <- find_cats_data(id)
  tryCatch({
    # Read time zone from tag guide
    tag_guide <- readxl::read_xlsx("/Volumes/COPYCATSdat/CATS/TAG GUIDE.xlsx", skip = 2)
    tag_tz <- sprintf("Etc/GMT%+d", -filter(tag_guide, ID == id)$`UTC         _`)
    prh <- readMat(data_paths$prh_path)
    tagonoff_prh <- dn_to_posix(range(prh$DN[prh$tagon == 1])) %>% 
      force_tz(tag_tz)
    # Tag_On seems to read correctly, but Tag_Off comes across as a serial
    tagonoff_guide <- filter(tag_guide, ID == id)$Tag_On %>% 
      force_tz(tag_tz)
    tagonoff_guide[2] <- filter(tag_guide, ID == id)$Tag_Off %>% 
      openxlsx::convertToDateTime(tz = tag_tz)
    tibble(tagon_prh = tagonoff_prh[1],
           tagoff_prh = tagonoff_prh[2],
           tagon_guide = tagonoff_guide[1],
           tagoff_guide = tagonoff_guide[2])
  }, error = function(e) tibble(tagon_prh = NA,
                                tagoff_prh = NA,
                                tagon_guide = NA,
                                tagoff_guide = NA))
}

# Deployment metadata. Starts with Paolo's lunge rate file. Filters down to
# good quality dives and the four rorqual species of interest. Removes sonar
# exposure deployments. Gets tag on/off times.
deployments_cats <- read_csv("data/lunge_rates_from_Paolo.csv") %>% 
  filter(lunge_quality %in% c("good", "good dives", "good_dives"), 
         sonar_exp == "none",
         species %in% c("ba", "bp", "bw", "mn")) %>% 
  group_by(species, ID) %>% 
  group_map(tag_times_cats)
deployments_dtag <- read_xlsx("/Volumes/COPYCATSdat/DTag data.xlsx", 
                              col_types = c("text", "date", "date", "numeric", 
                                            "numeric", "text", "numeric", 
                                            "numeric", "text", "text", 
                                            "text")) %>% 
  mutate(species = substr(ID, 1, 2),
         cee_start = ISOdatetime(year(`Tag on`), 
                                 month(`Tag on`), 
                                 day(`Tag on`), 
                                 ifelse(nchar(`Time CEE start`) == 4, 
                                        substr(`Time CEE start`, 1, 2), 
                                        substr(`Time CEE start`, 1, 1)), 
                                 ifelse(nchar(`Time CEE start`) == 4, 
                                        substr(`Time CEE start`, 3, 4), 
                                        substr(`Time CEE start`, 2, 3)), 
                                 0, 
                                 tz = "UTC"), 
         tagon_ceil = ceiling_date(`Tag on`, "hours"), 
         tagend = coalesce(cee_start, `Tag off`), 
         tagoff_floor = floor_date(tagend, "hours"), 
         data_hours = as.numeric(tagoff_floor - tagon_ceil, "hours"))

# Looks up the lunge times for all the deployments.
lunge_tbl_cats <- deployments_cats %>% 
  group_by_all %>% 
  group_map(lunge_times_cats)
lunge_tbl_dtag <- deployments_dtag %>% 
  group_by_all %>% 
  group_map(lunge_times_dtag)

lunge_tbl <- ungroup(lunge_tbl_cats) %>% 
  rbind(select(ungroup(lunge_tbl_dtag),
               species,
               ID,
               tagon = `Tag on`,
               tagoff = `Tag off`,
               lunge_dt))
 
# Final product. A tibble with columns for whale ID, species, tag on/off time, 
# hour of day, and feeding rate for that hour. Fills in missing hours.
# Drops leading and trailing incomplete hours. For example, a deployment from
# 12:35 to 14:20 would contribute one row: 13:00 - 14:00.
hourly_lunges <- lunge_tbl %>% 
  drop_na() %>%
  mutate(date = date(lunge_dt),
         hour = hour(lunge_dt),
         start = ceiling_date(tagon, "hours"),
         end = floor_date(tagoff, "hours") - hours(1)) %>% 
  filter(lunge_dt >= start,
         lunge_dt <= end) %>% 
  group_by_at(vars(-lunge_dt)) %>% 
  summarize(rf_h = n()) %>% 
  group_by_at(vars(-date, -hour, -rf_h)) %>% 
  group_map(function(data, key) {
    tibble(bucket = seq(key$start, key$end, "hours"),
           date = date(bucket),
           hour = hour(bucket)) %>% 
      left_join(data, by = c("date", "hour")) %>% 
      replace_na(list(rf_h = 0))
  }) %>% 
  ungroup %>% 
  mutate(binomial = recode(species, 
                           ba = "Balaenoptera acutorostrata",
                           bp = "Balaenoptera physalus",
                           bw = "Balaenoptera musculus",
                           mn = "Megaptera novaeangliae"))

# Box plots of hourly feeding rates
hourly_lunges %>% 
  group_by(species, binomial) %>% 
  group_walk(function(data, key) {
    boxplot <- data %>% 
      ggplot(aes(x = hour, y = rf_h, group = hour)) +
      geom_boxplot() +
      labs(y = "Lunges",
           title = key$binomial) +
      theme_classic() +
      theme(axis.title.x = element_blank())
    nplot <- data %>% 
      ggplot(aes(x = hour)) +
      geom_bar() +
      labs(x = "Hour of day",
           y = "Sample size") +
      theme_classic()
    boxplot + nplot + patchwork::plot_layout(ncol = 1, heights = c(4, 1))
    ggsave(sprintf("figs/FORMATT/%s.pdf", key$species), height = 6, width = 9)
  })
