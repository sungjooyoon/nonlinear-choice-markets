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


```r
# Prelec function
prelec_w <- function(q, gamma, delta) {
  exp(-delta * (-log(q))^gamma)
}

# Fitting Prelec function
fit_prelec <- function(data) {
  nls(actual_prob ~ prelec_w(implied_prob, gamma, delta), data = data,
      start = list(gamma = 1, delta = 1), algorithm = "port")
}

# Intervals as defined earlier in a structure
intervals <- c("first_implied", "last_15_to_12", "last_12_to_9", "last_9_to_6", "last_6_to_3", "last_3_to_0")

# Prepare data for fitting
fit_results <- list()

# Fit Prelec function for each interval
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
  prelec_fit <- fit_prelec(subset_data)
  gamma <- coef(prelec_fit)["gamma"]
  delta <- coef(prelec_fit)["delta"]
  fit_results[[interval]] <- list(interval = interval, gamma = gamma, delta = delta, fit = prelec_fit, data = subset_data)
}

# Extract fitted parameters for each interval and compile into a data frame
fit_parameters <- map_dfr(fit_results, function(res) {
  data.frame(
    interval = res$interval,
    gamma = res$gamma,
    delta = res$delta,
    observations = nrow(res$data),
    clusters = n_distinct(res$data$event_id),
    p_value = summary(res$fit)$coefficients["gamma", "Pr(>|t|)"]
  )
})

# Print the fitted parameters table
fit_parameters |>
  kable(col.names = c("Interval", "Gamma (Curvature)", "Delta (Elevation)", 
                      "Observations", "Clusters", "p-value (H0: Gamma = 1)"))
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


```r
# Gamma paramerer over time
fit_parameters$interval <- factor(fit_parameters$interval, levels = c("first_implied", "last_15_to_12", "last_12_to_9", "last_9_to_6", "last_6_to_3", "last_3_to_0"))
ggplot(fit_parameters, aes(x = interval, y = gamma)) +
  geom_line() +
  geom_point() +
  labs(title = "Gamma Parameter Over Time",
       x = "Time Interval",
       y = "Gamma") +
  theme_minimal()
```

```
## `geom_line()`: Each group consists of only one observation.
## ℹ Do you need to adjust the group aesthetic?
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

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
