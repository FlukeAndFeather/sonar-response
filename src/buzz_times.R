library(tidyverse)

# Odontocete feeding rate quantiles
## Deployment metadata
buzz_tbl <- read_csv("data/deployment_info.csv") %>% 
  # Remove CEEs
  filter(exposure == 0) %>% 
  mutate_at(vars(`start time UTC`, `end time UTC`), lubridate::mdy_hm) %>% 
  # Buzz times
  left_join(read_csv("data/all_buzzes.csv")) %>% 
  # Clean up species binomial name
  mutate(binomial = str_replace(species, "_", " "), 
         # Calculate hour of deployment
         buzz_dt = `start time UTC` + buzz_start,
         begin = `start time UTC` + (as.numeric(`end time UTC` - `start time UTC`, units = "secs") %% 3600),
         hour = floor(as.numeric(buzz_dt - begin, units = "hours"))) %>% 
    # Remove the first hour of Pp deployments
    filter(!(binomial == "Phocoena phocoena" & hour == 0)) %>% 
    select(binomial, deployment_ID, hour, `start time UTC`, `end time UTC`, begin) %>% 
    # Group buzzes by hour of deployment
    group_by_all %>% 
    # Count buzzes per hour to get hourly feeding rates
    summarize(rf_h = n()) %>% 
    # Divide orca feeding rates by six (see Danuta's email)
    mutate(rf_h = if_else(binomial == "Orcinus orca",
                          as.integer(ceiling(rf_h / 6)),
                          rf_h)) %>% 
    group_by_at(vars(-hour, -rf_h)) %>% 
    # Fill in gaps
    group_map(function(data, key) {
      dep_hrs <- floor(as.numeric(key$`end time UTC` - key$`start time UTC`, 
                                  units = "hours"))
      tibble(hour = seq_len(dep_hrs)) %>% 
        left_join(data, by = "hour") %>% 
        replace_na(list(rf_h = 0))
    })
  
  # Quantile functions and ECDF plots of feeding rates per species
buzz_rf <- buzz_tbl %>% 
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
    ggsave(sprintf("figs/rf_cdfs/%s.pdf", binom),
           plot,
           width = 9,
           height = 6)
    # Quantile function
    q_buzz <- function(p, ...) {
      quantile(data$rf_h, probs = p)
    }
    
    tibble(mean_rf = mean(data$rf_h),
           var_rf = var(data$rf_h),
           median_rf = median(data$rf_h),
           firstq_rf = q_buzz(0.25),
           thirdq_rf = q_buzz(0.75),
           q_rf_fun = list(q_buzz))
  }) %>% 
  ungroup

save(buzz_rf, file = "data/buzz_rf.RData")
