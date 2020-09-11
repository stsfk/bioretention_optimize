# Purpose:
#
# The script in this file generates random rainfalls that are used to drive SWMM simulation.
# 
# The rainfalls are considered as superposition of some "base rainfalls".
# The base rainfall has a single peak (triangular shape), and the duration and the peak intensity are random.
# The base rainfall can happen at any time step at random.

# The random rainfall generation process is controlled by several parameters:
# 1. n: the number of base rainfalls
# 2. mu: mean of peak rainfall intensity of a base rainfall, assuming normal distribution
# 3. sigma: standard deviation of peak rainfall of a base rainfall, assuming normal distribution
# 4. l: length of the output rainfall time series
# 5. max_duration, minimal_duration: lower and upper limits of duration of base rainfalls

# Libraries ---------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, caret, lubridate, RcppRoll, zeallot, xgboost)



# Parameters --------------------------------------------------------------


set.seed(10000)

n <- 35  # number of base rainfalls
mu <- 0.5
sigma <- 0.7
l <- 365*24*6 # 10-min resolution rainfall time series for 1 year
max_duration <- 84
minimal_duration <- 5

base_shape <- c(1,2,6,2,1)


# Processing --------------------------------------------------------------

for (iter in 1:10){
  # repeat 10 times
  
  rain_locations <- runif(n, min = max_duration, max = l - max_duration) %>% round # assuming dry periods at begining and end of time series
  rain_durations <- runif(n, min = minimal_duration, max = max_duration) %>% round
  rain_peaks <- rnorm(n, mu, sigma) %>% abs()
  rain_peaks[rain_peaks > 5] <- 5 # peak rainfall of a base rainfall < 5
  
  out <- rep(0, l)
  for (i in 1:n){
    
    rain_series <- approx(x = 1:5, base_shape, n = rain_durations[i])$y/6*rain_peaks[i]
    
    out[rain_locations[i]:(rain_locations[i] + rain_durations[i] - 1)] <- out[rain_locations[i]:(rain_locations[i] + rain_durations[i] - 1)]  +
      rain_series
  }
  
  out <- tibble(
    datetime = seq(ymd_hm("2007-01-01 00:00"), by = 600, length.out = length(out)),
    rain = out
  )
  
  # visualize rainfall time series
  ggplot(out, aes(datetime, rain)) +
    geom_line() +
    labs(y = "Rainfall intensity [Inch/hour]")
  
  # write 
  write.table(
    out %>%
      mutate(lines = paste("rain", year(datetime), month(datetime), day(datetime), hour(datetime), minute(datetime), rain)) %>%
      pull(lines),
    file = paste0("./data/rain_", iter,".csv"),
    row.names = F, 
    quote = F, 
    col.names = F
  )
}