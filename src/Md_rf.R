library(tidyverse)

Md_deployments <- tribble(~whale,   ~date,        ~tagdur_hr, ~fd, ~fddur_min, ~fdmaxdep_m, ~socdep_m, ~search_min, ~buzz_ct,
                          "MdH1",   "2008-05-16", 18.4,       7,   48,         911,         448,       28,          35,
                          "MdH1",   "2005-10-21", 4.1,        3,   50,         671,         475,       22,          24,
                          "MdH1",   "2003-10-11", 12.5,       5,   51,         616,         414,       26,          26,
                          "MdH15",  "2003-10-25", 2.6,        2,   47,         774,         426,       25,          23,
                          "MdH22",  "2008-10-15", 18.0,       7,   44,         710,         340,       23,          21,
                          "MdH22",  "2005-10-21", 2.8,        1,   47,         616,         520,       21,          18,
                          "MdH22",  "2004-10-13", 9.5,        4,   44,         1003,        473,       28,          32,
                          "MdH6",   "2008-05-15", 2.0,        2,   48,         781,         389,       24,          23,
                          "MdH6",   "2005-10-04", 6.9,        3,   57,         914,         518,       25,          29,
                          "MdH43",  "2005-10-12", 8.6,        4,   45,         833,         505,       25,          43,
                          "MdH74",  "2008-05-21", 1.6,        1,   47,         807,         419,       20,          11,
                          "MdHC1",  "2008-05-27", 6.2,        2,   58,         932,         461,       27,          34,
                          "MdHX33", "2010-05-26", 2.9,        1,   48,         925,         503,       22,          34,
                          "MdH86",  "2010-06-10", 15.3,       8,   41,         834,         353,       22,          18) %>% 
  mutate(date = lubridate::ymd(date),
         id = sprintf("%s%s", format(date, "%y%m%d"), whale))

expand_dive <- function(row) {
  # Convert durations to seconds
  tagdur_s <- row$tagdur_hr * 60 * 60
  fddur_s <- row$fddur_min * 60
  cycle_s <- tagdur_s / row$fd
  
  # 1s sampling interval
  t <- 0:tagdur_s
  
  # Assume constant descent rate from surface to max depth
  desc_rate <- row$fdmaxdep_m / fddur_s / 2
  # Append a surface interval after the dive
  divedepth_m <- c(seq(0, row$fdmaxdep_m, length.out = fddur_s / 2),
                   seq(row$fdmaxdep_m, 0, length.out = fddur_s / 2),
                   rep(0, times = cycle_s - fddur_s + 1))
  # Repeat the dives
  depth_m <- rep(divedepth_m, row$fd)
  # Distribute buzzes below "start of clicking" depth
  soc_start <- first(t[which(divedepth_m > row$socdep_m)])
  soc_end <- last(t[which(divedepth_m > row$socdep_m)])
  buzz_t <- map(length(divedepth_m) * 0:(row$fd - 1),
                ~ .x + seq(soc_start, soc_end, length.out = row$buzz_ct)) %>% 
    unlist
  buzz_m <- depth_m[buzz_t]
  list(metadata = row,
       profile = tibble(t, depth_m = depth_m[1:length(t)]),
       buzzes = tibble(buzz_t, buzz_m))
}

Md_buzzes <- Md_deployments %>% 
  group_by(id) %>% 
  group_split() %>% 
  map(expand_dive)

hourly_buzzes <- function(buzzes, tagdur_hr) {
  buzzes %>% 
    mutate(hour = ceiling(buzz_t / 3600)) %>% 
    group_by(hour) %>% 
    summarize(rf = n()) %>% 
    right_join(tibble(hour = 1:floor(tagdur_hr)),
               by = "hour") %>% 
    replace_na(list(rf = 0))
}

buzz_plot <- function(buzz) {
  caption <- with(buzz$metadata,
                  sprintf("tagdur_hr = %.1f, fd = %d, buzz_ct = %d", 
                          tagdur_hr, fd, buzz_ct))
  hourly <- hourly_buzzes(buzz$buzzes, buzz$metadata$tagdur_hr)
  
  ggplot() +
    geom_line(aes(t, depth_m),
              buzz$profile) +
    geom_point(aes(buzz_t, buzz_m),
               buzz$buzzes) +
    geom_vline(xintercept = seq(3600, floor(buzz$metadata$tagdur_hr) * 3600, by = 3600),
               linetype = "dashed") +
    geom_text(aes(x, y, label = label),
              tibble(x = seq(1800, floor(buzz$metadata$tagdur_hr) * 3600, by = 3600),
                     y = -25,
                     label = hourly$rf)) +
    scale_y_reverse() +
    labs(x = "Seconds since start",
         y = "Depth (m)",
         title = buzz$metadata$id,
         caption = caption) +
    theme_minimal()
}

Md_buzz_plots <- Md_buzzes %>% 
  map(buzz_plot)

Md_tbl <- Md_buzzes %>% 
  map_dfr(~ cbind(id = .x$metadata$id,
                  hourly_buzzes(.x$buzzes, .x$metadata$tagdur_hr)))

# Md rf ECDF
ggplot(Md_tbl, aes(x = rf)) +
  stat_ecdf() +
  theme_minimal()

# Md rf (made to match buzz_rf, see buzz_times.R)
q_buzz <- function(p, ...) {
  quantile(Md_tbl$rf, probs = p)
}
Md_buzz_rf <- Md_tbl %>% 
  mutate(binomial = "Mesoplodon densirostris") %>% 
  group_by(binomial) %>% 
  summarize(mean_rf = mean(rf),
            sd_rf = sd(rf),
            median_rf = median(rf),
            firstq_rf = q_buzz(0.25),
            thirdq_rf = q_buzz(0.75),
            q_rf_fun = list(q_buzz))

save(Md_buzz_rf, file = "data/Md_buzz_rf.RData")
