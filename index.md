---
title: "Nonlinear Prob. Weighting, draft"
author: "Sungjoo Yoon — Camerer Group @ Caltech"
date: "2024-07-06"
output:
  distill::distill_article:
---

# Setup
To play with the QJE data, I recommend installing the 'readstata13' package hosted on CRAN. Other packages like 'haven' are notoriously slow.
If your computer meets the requirements (64-bit Windows OS, sufficient RAM), I would also suggest manually increasing memory alloc.
While 'include=FALSE' is on so this is not visible, I'm using the 'dplyr', 'ggplot2', 'knitr', 'markdown', 'modelsummary', 'nls2', 'purrr', 'readstata13', and 'tidyverse' packages.


# Notes

- Goal is to see how choice and risk changes as priors stack up towards end of timed market
- Prelec one parameter: $$ w(p) = \exp(-(-\ln(p))^\alpha) $$
- Prelec two parameter: $$ w(p) = \exp(-\beta(-\ln(p))^\alpha) $$
- High-level if we define these functions in code and run NLS to fit data, we can see behavior of gamma curvature (actual vs. implied probs)
- Gamma curvature reflects risk attitude / how much weighting function distorts perceived vs. objective
- Prof. Camerer: if gamma is just >1, S-shaped where pi(p)<p. 
- This is opposite of the typical pattern where pi(p)>p for p<1/e and pi(p)<p for p>1/e
- "In Prelec one parameter form pi(1/e)=1/e is identity point around which pi(p) rotates"

# Data

Running this next chunk will reflect why I recommend memory increases.


```r
# Reading in the data (takes an unholy amount of system memory)
dataFull <- read.dta13("QJE_data_reconstructed/betfair/Data/dataCombined2.dta")
```


```r
# Reformatting and filtering to our downsample
dataFull <- dataFull |>
  filter(sports_id == "Basketball") |>
  mutate(
    dtactual_off = as.POSIXct(dtactual_off, format="%Y-%m-%d %H:%M:%S"),
    first_taken = as.POSIXct(first_taken, format="%Y-%m-%d %H:%M:%S")
  )

# Seed for repro — shoutout Pasadena
set.seed(91125)
dataFull <- dataFull |>
  sample_n(100000)

# Function to process data chunks and mutate timestamps
process_chunk <- function(data_chunk) {
  data_chunk <- data_chunk |>
    mutate(
      event_duration = as.numeric(difftime(first_taken, dtactual_off, units = "secs"))
    ) |>
    group_by(event_id) |>
    mutate(
      event_total_duration = max(event_duration, na.rm = TRUE)
    ) |>
    ungroup()
  return(data_chunk)
}
chunk_size <- 10000 
data_chunks <- split(dataFull, (seq(nrow(dataFull)) - 1) %/% chunk_size)
processed_chunks <- lapply(data_chunks, process_chunk)
dataFull_processed <- bind_rows(processed_chunks)
```

Note: pre-game odds are not available in this dataset, exclusively in-play. So instead, I took the earliest in-play bet available. Please let me know if you want me to change the approach. Otherwise, the following code also allows us to segment the last fifteen minutes of the event.


```r
# Calculating the earliest in-play bet available, which also serves as first implied odds
earliest_bets <- dataFull_processed |>
  group_by(event_id) |>
  filter(first_taken == min(first_taken)) |>
  select(event_id, first_implied = perc) |>
  ungroup()

# Merge this information back into the main dataframe
dataFull_processed <- dataFull_processed |>
  left_join(earliest_bets, by = "event_id")
```

```
## Warning in left_join(dataFull_processed, earliest_bets, by = "event_id"): Detected an unexpected many-to-many relationship between `x` and `y`.
## ℹ Row 61 of `x` matches multiple rows in `y`.
## ℹ Row 34 of `y` matches multiple rows in `x`.
## ℹ If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.
```

```r
# Define the late segment intervals and create time_interval variable
dataFull_processed <- dataFull_processed |>
  filter(event_duration >= (event_total_duration - 15*60)) |>
  mutate(
    time_interval = case_when(
      event_duration >= (event_total_duration - 15*60) & event_duration < (event_total_duration - 12*60) ~ "last_15_to_12",
      event_duration >= (event_total_duration - 12*60) & event_duration < (event_total_duration - 9*60) ~ "last_12_to_9",
      event_duration >= (event_total_duration - 9*60) & event_duration < (event_total_duration - 6*60) ~ "last_9_to_6",
      event_duration >= (event_total_duration - 6*60) & event_duration < (event_total_duration - 3*60) ~ "last_6_to_3",
      event_duration >= (event_total_duration - 3*60) ~ "last_3_to_0"
    )
  )
```

# Analysis

Prelec function and fitting

Draft 1 Not sure if this is right

```r
# Function for calibration curve
calculate_calibration_curve <- function(data, bins = 10) {
  data <- data |>
    mutate(bin = cut(implied_prob, breaks = bins, labels = FALSE)) |>
    group_by(bin) |>
    summarize(
      mean_pred_prob = mean(implied_prob),
      mean_actual_prob = mean(actual_prob)
    )
  return(data)
}

# Prelec 1 parameter function
prelec_w_one_param <- function(q, gamma) {
  exp(-(-log(q))^gamma)
}

# Prelec 2 parameter function
prelec_w_two_param <- function(q, gamma, delta) {
  exp(-delta * (-log(q))^gamma)
}

# Fitting Prelec 1 parameter function to calibration curve
fit_prelec_one_param_calibration <- function(data) {
  nls(mean_actual_prob ~ prelec_w_one_param(mean_pred_prob, gamma), data = data,
      start = list(gamma = 1), algorithm = "port")
}

# Fitting Prelec 2 parameter function to calibration curve
fit_prelec_two_param_calibration <- function(data) {
  nls(mean_actual_prob ~ prelec_w_two_param(mean_pred_prob, gamma, delta), data = data,
      start = list(gamma = 1, delta = 1), algorithm = "port")
}

# Function to calculate calibration curves
calculate_calibration_curve <- function(data, bins = 10) {
  data <- data |>
    mutate(bin = cut(implied_prob, breaks = bins, labels = FALSE)) |>
    group_by(bin) |>
    summarize(
      mean_pred_prob = mean(implied_prob),
      mean_actual_prob = mean(actual_prob)
    )
  return(data)
}

# Intervals as defined earlier in a structure
intervals <- c("first_implied", "last_15_to_12", "last_12_to_9", "last_9_to_6", "last_6_to_3", "last_3_to_0")

# Prepare data for fitting with calibration curves
fit_results_one_param_calibration <- list()
fit_results_two_param_calibration <- list()

# Fit Prelec function for each interval (including first implied) w/ cal.
for (interval in intervals) {
  if (interval == "first_implied") {
    subset_data <- dataFull_processed |>
      group_by(event_id) |>
      filter(first_taken == min(first_taken)) |>
      mutate(implied_prob = perc, actual_prob = win_flag)
  } else {
    subset_data <- dataFull_processed |>
      filter(time_interval == interval) |>
      mutate(implied_prob = perc, actual_prob = win_flag)
  }

  # Calculate calibration curve
  calibration_data <- calculate_calibration_curve(subset_data)

  # Fit Prelec one-parameter model to calibration curve
  prelec_fit_one_param_calibration <- fit_prelec_one_param_calibration(calibration_data)
  gamma_one_param_calibration <- coef(prelec_fit_one_param_calibration)["gamma"]
  fit_results_one_param_calibration[[interval]] <- list(interval = interval, gamma = gamma_one_param_calibration, fit = prelec_fit_one_param_calibration, data = calibration_data)

  # Fit Prelec two-parameter model to calibration curve
  prelec_fit_two_param_calibration <- fit_prelec_two_param_calibration(calibration_data)
  gamma_two_param_calibration <- coef(prelec_fit_two_param_calibration)["gamma"]
  delta_two_param_calibration <- coef(prelec_fit_two_param_calibration)["delta"]
  fit_results_two_param_calibration[[interval]] <- list(interval = interval, gamma = gamma_two_param_calibration, delta = delta_two_param_calibration, fit = prelec_fit_two_param_calibration, data = calibration_data)
}

# Fitted parameters to new df structure for 1p model using calibration curves
fit_parameters_one_param_calibration <- map_dfr(fit_results_one_param_calibration, function(res) {
  data.frame(
    interval = res$interval,
    gamma = res$gamma,
    observations = nrow(res$data),
    clusters = n_distinct(res$data$event_id),
    p_value = summary(res$fit)$coefficients["gamma", "Pr(>|t|)"]
  )
})
```

```
## Warning: Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
```

```r
# Fitted parameters to new df structure for 2p model using calibration curves
fit_parameters_two_param_calibration <- map_dfr(fit_results_two_param_calibration, function(res) {
  data.frame(
    interval = res$interval,
    gamma = res$gamma,
    delta = res$delta,
    observations = nrow(res$data),
    clusters = n_distinct(res$data$event_id),
    p_value = summary(res$fit)$coefficients["gamma", "Pr(>|t|)"]
  )
})
```

```
## Warning: Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
## Unknown or uninitialised column: `event_id`.
```

```r
# Tables
fit_parameters_one_param_calibration |>
  kable(col.names = c("Interval", "Gamma (Curvature)", "Observations", "Clusters", "p-value (H0: Gamma = 1)"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Interval </th>
   <th style="text-align:right;"> Gamma (Curvature) </th>
   <th style="text-align:right;"> Observations </th>
   <th style="text-align:right;"> Clusters </th>
   <th style="text-align:right;"> p-value (H0: Gamma = 1) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> gamma...1 </td>
   <td style="text-align:left;"> first_implied </td>
   <td style="text-align:right;"> 1.5876639 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0084363 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...2 </td>
   <td style="text-align:left;"> last_15_to_12 </td>
   <td style="text-align:right;"> 0.9601368 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000564 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...3 </td>
   <td style="text-align:left;"> last_12_to_9 </td>
   <td style="text-align:right;"> 0.8087471 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0001752 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...4 </td>
   <td style="text-align:left;"> last_9_to_6 </td>
   <td style="text-align:right;"> 1.1587038 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000041 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...5 </td>
   <td style="text-align:left;"> last_6_to_3 </td>
   <td style="text-align:right;"> 1.2266805 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000002 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...6 </td>
   <td style="text-align:left;"> last_3_to_0 </td>
   <td style="text-align:right;"> 1.0558439 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
</tbody>
</table>

```r
fit_parameters_two_param_calibration |>
  kable(col.names = c("Interval", "Gamma (Curvature)", "Delta (Elevation)", "Observations", "Clusters", "p-value (H0: Gamma = 1)"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Interval </th>
   <th style="text-align:right;"> Gamma (Curvature) </th>
   <th style="text-align:right;"> Delta (Elevation) </th>
   <th style="text-align:right;"> Observations </th>
   <th style="text-align:right;"> Clusters </th>
   <th style="text-align:right;"> p-value (H0: Gamma = 1) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> gamma...1 </td>
   <td style="text-align:left;"> first_implied </td>
   <td style="text-align:right;"> 1.5862352 </td>
   <td style="text-align:right;"> 0.9988483 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.2179916 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...2 </td>
   <td style="text-align:left;"> last_15_to_12 </td>
   <td style="text-align:right;"> 0.9142696 </td>
   <td style="text-align:right;"> 0.8743932 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000871 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...3 </td>
   <td style="text-align:left;"> last_12_to_9 </td>
   <td style="text-align:right;"> 0.7076741 </td>
   <td style="text-align:right;"> 0.7772526 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000133 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...4 </td>
   <td style="text-align:left;"> last_9_to_6 </td>
   <td style="text-align:right;"> 1.1023219 </td>
   <td style="text-align:right;"> 0.9026284 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000073 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...5 </td>
   <td style="text-align:left;"> last_6_to_3 </td>
   <td style="text-align:right;"> 1.1805037 </td>
   <td style="text-align:right;"> 0.9196380 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000003 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...6 </td>
   <td style="text-align:left;"> last_3_to_0 </td>
   <td style="text-align:right;"> 1.1026350 </td>
   <td style="text-align:right;"> 1.0938470 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
</tbody>
</table>

```r
# Visualization for gamma parameter in two-parameter model
fit_parameters_one_param_calibration$interval <- factor(fit_parameters_one_param_calibration$interval, levels = intervals)
fit_parameters_two_param_calibration$interval <- factor(fit_parameters_two_param_calibration$interval, levels = intervals)

ggplot(fit_parameters_one_param_calibration, aes(x = interval, y = gamma)) +
  geom_line(group = 1) +
  geom_point() +
  labs(title = "Gamma Parameter for One-Parameter Model",
       x = "Interval",
       y = "Gamma (Curvature)")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
ggplot(fit_parameters_two_param_calibration, aes(x = interval, y = gamma)) +
  geom_line(group = 1) +
  geom_point() +
  labs(title = "Gamma Parameter for Two-Parameter Model",
       x = "Interval",
       y = "Gamma (Curvature)")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)

Draft 0.5: just using regression to fit the data

```r
# Prelec 1 parameter function
prelec_w_one_param <- function(q, gamma) {
  exp(-(-log(q))^gamma)
}

# Fitting Prelec 1 parameter function
fit_prelec_one_param <- function(data) {
  nls(actual_prob ~ prelec_w_one_param(implied_prob, gamma), data = data,
      start = list(gamma = 1), algorithm = "port")
}

# Prelec 2
prelec_w_two_param <- function(q, gamma, delta) {
  exp(-delta * (-log(q))^gamma)
}

# Fitting Prelec 2 parameter function
fit_prelec_two_param <- function(data) {
  nls(actual_prob ~ prelec_w_two_param(implied_prob, gamma, delta), data = data,
      start = list(gamma = 1, delta = 1), algorithm = "port")
}

# Intervals as defined earlier in a structure
intervals <- c("first_implied", "last_15_to_12", "last_12_to_9", "last_9_to_6", "last_6_to_3", "last_3_to_0")

# Prepare data for fitting
fit_results_one_param <- list()
fit_results_two_param <- list()

# Fit Prelec function for each interval (including first implied)
for (interval in intervals) {
  if (interval == "first_implied") {
    subset_data <- dataFull_processed |>
      group_by(event_id) |>
      filter(first_taken == min(first_taken)) |>
      mutate(implied_prob = perc, actual_prob = win_flag)
  } else {
    subset_data <- dataFull_processed |>
      filter(time_interval == interval) |>
      mutate(implied_prob = perc, actual_prob = win_flag)
  }

  # Fit Prelec one-parameter model
  prelec_fit_one_param <- fit_prelec_one_param(subset_data)
  gamma_one_param <- coef(prelec_fit_one_param)["gamma"]
  fit_results_one_param[[interval]] <- list(interval = interval, gamma = gamma_one_param, fit = prelec_fit_one_param, data = subset_data)

  # Fit Prelec two-parameter model
  prelec_fit_two_param <- fit_prelec_two_param(subset_data)
  gamma_two_param <- coef(prelec_fit_two_param)["gamma"]
  delta_two_param <- coef(prelec_fit_two_param)["delta"]
  fit_results_two_param[[interval]] <- list(interval = interval, gamma = gamma_two_param, delta = delta_two_param, fit = prelec_fit_two_param, data = subset_data)
}

# Fitted parameters to new df structure for one-parameter model
fit_parameters_one_param <- map_dfr(fit_results_one_param, function(res) {
  data.frame(
    interval = res$interval,
    gamma = res$gamma,
    observations = nrow(res$data),
    clusters = n_distinct(res$data$event_id),
    p_value = summary(res$fit)$coefficients["gamma", "Pr(>|t|)"]
  )
})

# Fitted parameters to new df structure for two-parameter model
fit_parameters_two_param <- map_dfr(fit_results_two_param, function(res) {
  data.frame(
    interval = res$interval,
    gamma = res$gamma,
    delta = res$delta,
    observations = nrow(res$data),
    clusters = n_distinct(res$data$event_id),
    p_value = summary(res$fit)$coefficients["gamma", "Pr(>|t|)"]
  )
})

# Tables
fit_parameters_one_param |>
  kable(col.names = c("Interval", "Gamma (Curvature)", "Observations", "Clusters", "p-value (H0: Gamma = 1)"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Interval </th>
   <th style="text-align:right;"> Gamma (Curvature) </th>
   <th style="text-align:right;"> Observations </th>
   <th style="text-align:right;"> Clusters </th>
   <th style="text-align:right;"> p-value (H0: Gamma = 1) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> gamma...1 </td>
   <td style="text-align:left;"> first_implied </td>
   <td style="text-align:right;"> 0.9883489 </td>
   <td style="text-align:right;"> 34991 </td>
   <td style="text-align:right;"> 33612 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...2 </td>
   <td style="text-align:left;"> last_15_to_12 </td>
   <td style="text-align:right;"> 0.9707036 </td>
   <td style="text-align:right;"> 683 </td>
   <td style="text-align:right;"> 649 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...3 </td>
   <td style="text-align:left;"> last_12_to_9 </td>
   <td style="text-align:right;"> 0.9544373 </td>
   <td style="text-align:right;"> 774 </td>
   <td style="text-align:right;"> 732 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...4 </td>
   <td style="text-align:left;"> last_9_to_6 </td>
   <td style="text-align:right;"> 1.2542790 </td>
   <td style="text-align:right;"> 871 </td>
   <td style="text-align:right;"> 823 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...5 </td>
   <td style="text-align:left;"> last_6_to_3 </td>
   <td style="text-align:right;"> 1.2874838 </td>
   <td style="text-align:right;"> 1050 </td>
   <td style="text-align:right;"> 980 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...6 </td>
   <td style="text-align:left;"> last_3_to_0 </td>
   <td style="text-align:right;"> 1.0487914 </td>
   <td style="text-align:right;"> 86092 </td>
   <td style="text-align:right;"> 33612 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

```r
fit_parameters_two_param |>
  kable(col.names = c("Interval", "Gamma (Curvature)", "Delta (Elevation)", "Observations", "Clusters", "p-value (H0: Gamma = 1)"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Interval </th>
   <th style="text-align:right;"> Gamma (Curvature) </th>
   <th style="text-align:right;"> Delta (Elevation) </th>
   <th style="text-align:right;"> Observations </th>
   <th style="text-align:right;"> Clusters </th>
   <th style="text-align:right;"> p-value (H0: Gamma = 1) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> gamma...1 </td>
   <td style="text-align:left;"> first_implied </td>
   <td style="text-align:right;"> 1.0229094 </td>
   <td style="text-align:right;"> 1.0697222 </td>
   <td style="text-align:right;"> 34991 </td>
   <td style="text-align:right;"> 33612 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...2 </td>
   <td style="text-align:left;"> last_15_to_12 </td>
   <td style="text-align:right;"> 0.8680689 </td>
   <td style="text-align:right;"> 0.8616903 </td>
   <td style="text-align:right;"> 683 </td>
   <td style="text-align:right;"> 649 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...3 </td>
   <td style="text-align:left;"> last_12_to_9 </td>
   <td style="text-align:right;"> 0.7743088 </td>
   <td style="text-align:right;"> 0.8037398 </td>
   <td style="text-align:right;"> 774 </td>
   <td style="text-align:right;"> 732 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...4 </td>
   <td style="text-align:left;"> last_9_to_6 </td>
   <td style="text-align:right;"> 1.1681824 </td>
   <td style="text-align:right;"> 0.9164937 </td>
   <td style="text-align:right;"> 871 </td>
   <td style="text-align:right;"> 823 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...5 </td>
   <td style="text-align:left;"> last_6_to_3 </td>
   <td style="text-align:right;"> 1.2196140 </td>
   <td style="text-align:right;"> 0.9331590 </td>
   <td style="text-align:right;"> 1050 </td>
   <td style="text-align:right;"> 980 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gamma...6 </td>
   <td style="text-align:left;"> last_3_to_0 </td>
   <td style="text-align:right;"> 1.1021232 </td>
   <td style="text-align:right;"> 1.0962072 </td>
   <td style="text-align:right;"> 86092 </td>
   <td style="text-align:right;"> 33612 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

Viz 0.5


```r
fit_parameters_one_param$interval <- factor(fit_parameters_one_param$interval, levels = intervals)
fit_parameters_two_param$interval <- factor(fit_parameters_two_param$interval, levels = intervals)

plot_one_param <- ggplot(fit_parameters_one_param, aes(x = interval, y = gamma)) +
  geom_point(size = 3) +
  geom_line(group = 1) +
  labs(title = "Prelec One-Parameter Model: Gamma Across Intervals",
       x = "Interval",
       y = "Gamma (Curvature)")

# Plot for two-parameter model
plot_two_param <- ggplot(fit_parameters_two_param, aes(x = interval, y = gamma)) +
  geom_point(size = 3) +
  geom_line(group = 1) +
  labs(title = "Prelec Two-Parameter Model: Gamma Across Intervals",
       x = "Interval",
       y = "Gamma (Curvature)")

plot_one_param
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

```r
plot_two_param
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png)

```r
# Viz on number of bets / sizes of trade
ggplot(dataFull_processed, aes(x = number_bets)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Number of Bets",
       x = "Number of Bets",
       y = "Frequency")
```

```
## Warning: Removed 2595 rows containing non-finite values (`stat_bin()`).
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-3.png)

Easy way to see the data at every step without bricking your computer, just run this chunk


```r
first_10000_rows <- head(dataFull_processed, 10000)
view(first_10000_rows)
```

GitHub Pages index.html generator script


```r
knit("index.Rmd", output = "index.md")
```

```
## 
## 
## processing file: index.Rmd
```

```
## Error in parse_block(g[-1], g[1], params.src, markdown_mode): Duplicate chunk label 'setup', which has been used for the chunk:
## knitr::opts_chunk$set(echo = TRUE)
## library(dplyr)
## library(ggplot2)
## library(knitr)
## library(markdown)
## library(modelsummary)
## library(nls2)
## library(purrr)
## library(readstata13)
## library(tidyverse)
```

```r
markdownToHTML("index.md", "index.html")
```

