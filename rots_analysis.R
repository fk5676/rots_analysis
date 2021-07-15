# required libraries
library(knitr)
library(tidyverse)
library(broom)

# summary statistics function for "uitloging" and "gidsparameter"
# values (values are rounded to two significant digits similar
# to the memo data)
summary_stats <- function(x){
    x %>%
    summarise_at(
    c("uitloging",
      "gidsparameter"),
       list(mean = mean,
            median = median,
            min = min,
            max = max), na.rm = TRUE) %>%
  round(digits = 2) %>%
  pivot_longer(cols = 1:ncol(.)) %>%
  separate(name, c("name", "statistic")) %>%
  arrange(name) %>%
  pivot_wider(id_cols = c(name, statistic)) %>%
  relocate(gidsparameter, .after = uitloging)
}

df <- readr::read_csv("data/data_rots.csv") %>%
  as_tibble()

# tabulate the data
knitr::kable(df,
             caption = "",
             col.names = c("Schudproef / kolomproef",
                           "monster nummer",
                           "Uitloging ($\\mu$g/l)",
                           "Concentratie gidsparameter ($\\mu$g/kg ds)"))

# simple stats
st <- df %>%
  summary_stats()

knitr::kable(st,
             col.names = c("Summary statistic",
                           "Uitloging ($\\mu$g/l)",
                           "Concentratie gidsparameter ($\\mu$g/kg ds)"))

# simple stats
st_limited <- df %>%
  filter(uitloging < 10) %>%
  summary_stats()

# read in rots parameters
# read in the data
parameters_rots <- readr::read_csv("data/parameters_rots.csv") %>%
  as_tibble()

# adds a column identifying the set as "complete"
df <- df %>%
  mutate(
    set = "full" 
  )

# reformat data for speedy evaluation
df <- df %>%
  filter(
    uitloging < 10
  ) %>%
  mutate(
    set = "limited"
  ) %>% 
  bind_rows(df) # bind output to the original input

# model data grouped by proef and set (limited or full dataset)
power <- df %>%
  group_by(proef, set) %>%
  do({
    data <- .
    fit <- lm(
      log(gidsparameter) ~ log(uitloging),
      data = data)
    tidy <- tidy(fit)
    tidy$r.squared <- glance(fit)$r.squared
    tidy
  })

# finally provide the combined results
exp <- df %>%
  group_by(proef, set) %>%
  do({
    data <- .
    fit <- lm(
      log(gidsparameter) ~ uitloging,
      data = data)
    tidy <- tidy(fit)
    tidy$r.squared <- glance(fit)$r.squared
    tidy
  })

# only group by set, to get the combined
# KP SP results
power_combined <- df %>%
  group_by(set) %>%
  do({
    data <- .
    fit <- lm(
      log(gidsparameter) ~ log(uitloging),
      data = data)
    tidy <- tidy(fit)
    tidy$r.squared <- glance(fit)$r.squared
    tidy
  })

# finally provide the combined results
exp_combined <- df %>%
  group_by(set) %>%
  do({
    data <- .
    fit <- lm(
      log(gidsparameter) ~ uitloging,
      data = data)
    tidy <- tidy(fit)
    tidy$r.squared <- glance(fit)$r.squared
    tidy
  })

# combine power and exponential results
power <- bind_rows(power, power_combined)
exp <- bind_rows(exp, exp_combined)

# replace terms for clarity
# and back-convert parameters from their linear form
power <- power %>%
  mutate(
    term = ifelse(grepl("uitloging",term),"b","a"),
    estimate = ifelse(term == "a", exp(estimate), estimate)
  ) %>%
  select(proef, set, term, estimate, r.squared) %>%
  pivot_wider(names_from = c(term),
              values_from = c(estimate, r.squared)) %>%
  mutate(
    func = "power"
  ) %>%
  rename(
    'r_squared' = 'r.squared_a',
    'a' = 'estimate_a',
    'b' = 'estimate_b'
  ) %>%
  select(-r.squared_b)

exp <- exp %>%
  mutate(
    term = ifelse(grepl("uitloging",term),"b","a"),
    estimate = ifelse(term == "a", exp(estimate), estimate)
  ) %>%
  select(proef, set, term, estimate, r.squared) %>%
  pivot_wider(names_from = c(term),
              values_from = c(estimate,r.squared)) %>%
  mutate(
    func = "exp"
  ) %>%
  rename(
    'r_squared' = 'r.squared_a',
    'a' = 'estimate_a',
    'b' = 'estimate_b'
  ) %>%
  select(-r.squared_b)

# bind the final parameters in a nice table
# clean up for reporting mainly relabelling
# combined "proef" results
parameters <- bind_rows(power, exp) %>%
  mutate(
    proef = ifelse(is.na(proef),"KP + SP", proef),
    a = round(a, 4),
    b = round(b, 4),
    r_squared = round(r_squared, 2)
  ) %>%
  ungroup() %>%
  arrange(set)

# combine both parameter sets for a comparison
parameters_diff <- dplyr::left_join(
  parameters,
  parameters_rots,
  by = c("proef", "set", "func"))

# take the difference of parameters and R2 values
parameters_diff <- parameters_diff %>%
  mutate(
    a_diff = a.x - a.y,
    b_diff = b.x - b.y,
    r_squared_diff = r_squared.x - r_squared.y
  ) %>%
  select(-c(ends_with(".y"))) %>%
  select(-c(ends_with(".x"))) %>%
  relocate(func, .after = r_squared_diff)

# set various standards
# on both parameter sets
parameters$uitloog_conc <- 4.5
parameters$epa_2016_conc <- 0.63
parameters$grondwater_conc <- 0.12
parameters$efsa_2020_conc <- 0.02

parameters_rots$uitloog_conc <- 4.5
parameters_rots$epa_2016_conc <- 0.63
parameters_rots$grondwater_conc <- 0.12
parameters_rots$efsa_2020_conc <- 0.02

# apply the models with fitted parameters
# to these standards, return in place
# (overwriting standard values)
parameters <- parameters %>%
  mutate(
    across(ends_with("_conc"),
    ~ round(
      ifelse(func == "power",
             a * .x ^ b,
             a * exp(.x * b))
      )
    )
  )

parameters_rots <- parameters_rots %>%
  mutate(
    across(ends_with("_conc"),
           ~ round(
             ifelse(func == "power",
                    a * .x ^ b,
                    a * exp(.x * b))
           )
    )
  )
