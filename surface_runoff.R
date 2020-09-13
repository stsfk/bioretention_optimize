## runoff reduction calculation

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, xts, zoo, lubridate)

# Read surface runoff time series -----------------------------------------
org_wd <- getwd()

x_gis_aggre <- list() # x_gis_aggre[[j]][[i]], j is year, i is different BC %
x_ud_aggre <- list() # x_ud[[j]][[i]]
x_org_aggre <- list() # orginal runoff without BC cover
x_runoff_aggre <- list()
for (j in 1:10){
  file_location <- str_c(org_wd, "/", "/results/year",j,"/",collapse = "")
  files_surf <- list.files(pattern='surf.*\\.txt', recursive=TRUE, path = file_location)
  files_surf <- paste("surf", c(1:(length(files_surf)-1)),".txt",sep="")
  files_ud <- paste("ud", c(1:(length(files_surf))),".txt",sep="")
  x_gis <- list()
  x_uds <- list()
  x_runoffs <- list()
  for (i in 1:81) {
    fid <- i
    file_name_surf <- paste0(file_location, files_surf[fid])
    file_name_ud <- paste0(file_location, files_ud[fid]) 
    x_gis[[i]] <- unname(unlist(read.table(file_name_surf, header = F)))/0.03 # convert from m3 per 1 min per 5000 m2 to L/s/ha
    x_uds[[i]] <- unname(unlist(read.table(file_name_ud, header = F)))/0.03
    x_runoffs[[i]] <- x_gis[[i]] + x_uds[[i]]
  }
  
  org_name <- paste0(file_location,"surf0.txt")
  x_org_aggre[[j]] <- unname(unlist(read.table(org_name)))/0.03
  
  x_gis_aggre[[j]] <- x_gis
  x_ud_aggre[[j]] <- x_uds
  x_runoff_aggre[[j]] <- x_runoffs
}

rm(list=setdiff(ls(), c("x_runoff_aggre","x_org_aggre")))  # clean

# Identify storm starts and ends ------------------------------------------
org_wd <- getwd()

read_rain <- function(file){
  temp <- readLines(file)
  
  rain_depths <- as.numeric(str_sub(temp, 18, -1))
  rain_depths<- rain_depths/12*25.4 # convert from inch/h to mm per 5 min
  
  return(c(rain_depths))
}

sta_end <- function(x, par.inter) {
  # This function is used to find start and end time of a storm event
  # Input:
  #   x is time series
  #   par.inter is the dry spell interval
  # Output:
  #   sta_end_data is a list containing the start end information
  loc.r <- which (x != 0)
  le <- length(loc.r)
  
  no.r <- rep(0, le)
  dur.r <- rep (0, le)
  
  j <- 1
  
  for (i in 2:le) {
    temp <- loc.r[i] - loc.r[i - 1]
    no.r[j] <- no.r[j] + 1
    if (temp > par.inter) {
      if (no.r[j] == 1) {
        dur.r[j] <- 1
      } else {
        dur.r[j] <- loc.r[i - 1] - loc.r[i - no.r[j]] + 1
      }
      j <- j + 1
    }
    
    if (i == le) {
      no.r[j] <- no.r[j] + 1
      dur.r[j] <-
        loc.r[length(loc.r)] - loc.r[length(loc.r) - no.r[j] + 1] + 1
    }
  }
  
  temp <- length(which(no.r != 0))
  no.r <- no.r[1:temp]
  dur.r <- dur.r[1:temp]
  
  rm(i, j, le, temp)
  
  m <- matrix(0, ncol = length(dur.r), nrow = max(dur.r))
  a <- ncol(m) # the number of rainfall event
  
  j = 1
  
  locmax <- rep(0, a)
  
  sta_end_data <- data.frame(sat = rep(0, a),
                             end = rep(0, a))
  for (i in 1:a) {
    temp <- loc.r[j] + dur.r[i] - 1
    m[1:dur.r[i], i] <- x[loc.r[j]:temp]
    sta_end_data[i, ] <- c(loc.r[j], temp)
    
    tempmax <- which.max(m[, i])
    locmax [i] <- loc.r[j] + tempmax - 1
    
    j <- j + no.r[i]
  }
  
  sta_end_data
}

setwd(paste0(org_wd,"/data"))
files <- list.files(pattern = '^rain.*\\.csv')
rains_list <- vector("list", length(files))

for (i in seq_along(files)){
  file <- paste0("rain_",i, ".csv")
  if (i != 1){
    rains_list[[i]] <- read_rain(file)    
  } else {
    rains_list[[i]] <- c(0,read_rain(file)) 
  }
}

inter <- as.list((1:4)*96)
sta_ends <- list() # sta_ends stores the start and end index for different par.inter
for (i in seq_along(rains_list)){
  rain <- c(rains_list[[i]])
  sta_ends[[i]] <- lapply(inter, sta_end, x = rain)
}

setwd(org_wd)

# Identify storm starts and ends 2nd part----------------------------------------------------
counter = 0

x_intervals <- list() # x_intervals is different storm interval
x_intervals_aggre <- list() # aggregated x_intervals
runoff_series_le <- length(x_org_aggre[[10]])

for (j in 1:10) {
  for (i in 1:length(inter)) {
    x_interval <- (sta_ends[[j]][[i]] - 1) * 5 + 1# every 5 mins, times 5 to convert to every minute, plus 1 is for correct indexing
    x_interval <- x_interval[, 1] # the start time
    x_intervals[[i]] <- c(0, x_interval, runoff_series_le)
  }
  x_intervals_aggre[[j]] <- x_intervals
}

# Discharge rate calculation --------------------------------------------------------

# Discharge rate calculation --------------------------------------------------------
find_peak_volume <- function(x, breaks){
  # This function is used to find runoff peak and volume for each individual storms
  labels <- 1:(length(breaks) - 1)
  
  x <- data.frame(x = x, id = seq_along(x)) %>%
    mutate(interval = cut(id, breaks=breaks, include.lowest = TRUE, right = F, labels = labels)) %>%
    group_by(interval) %>%
    summarise(volume = sum(x)/1000*60, # convert to m3 per ha
              peak = max(x)) 
  return(x)
}

find_peak_volumes <- function(x_aggre, x_intervals_aggre){
  peak_volumes <- vector("list", 10)
  for (i in 1:10){
    x <- x_aggre[[i]] # aggregated time series for 10 years
    peak_volume <- vector("list", 4) # 4 interval
    for (j in 1:4) {
      breaks <- x_intervals_aggre[[i]][[j]] # aggregated interval inds for different IETD
      temp <- find_peak_volume(x, breaks)
      temp$interval <- as.numeric(temp$interval)
      temp$year <- i
      temp$ietd <- j
      peak_volume[[j]] <- temp
    }
    peak_volume <- bind_rows(peak_volume)
    peak_volumes[[i]] <- peak_volume
  }
  temp <- bind_rows(peak_volumes)
  return(temp)
}

reduc_percentage <- function (x, y){
  (y - x)/y * 100
}

org_peak_volume <- 
  find_peak_volumes(x_aggre = x_org_aggre, x_intervals_aggre = x_intervals_aggre) # original peak volume

BR_peak_volume <- vector("list", 81) # Bioretention cases
for (i in 1:81){
  x_aggre <- vector("list", 10)
  for (j in 1:10){
    x_aggre[[j]] <- x_runoff_aggre[[j]][[i]]
  }
  
  df <- find_peak_volumes(x_aggre = x_aggre, x_intervals_aggre = x_intervals_aggre)
  
  df$gi <- i * 0.1 + 1.9
  
  BR_peak_volume[[i]] <- df
}

BR_peak_volumes <- bind_rows(BR_peak_volume)
BR_peak_volumes <- BR_peak_volumes %>%
  left_join(org_peak_volume, by = c("year","ietd","interval")) %>%
  select(ietd, year, interval, gi, everything()) %>%
  rename(volume = volume.x,
         peak = peak.x,
         org_volume = volume.y,
         org_peak = peak.y) %>%
  filter(org_volume != 0) %>%
  mutate(volume_perc = map2_dbl(volume, org_volume, reduc_percentage),
         peak_perc = map2_dbl(peak, org_peak, reduc_percentage),
         volume_perc = ifelse(volume_perc < 0, 0, volume_perc),
         peak_perc = ifelse(peak_perc < 0, 0, peak_perc)) 

# Figure 7 ----------------------------------------------------------------

data_plot <- BR_peak_volumes %>%
  filter(ietd == 2) %>%
  group_by(ietd, gi) %>%
  summarise(volume = sum(volume),
            org_volume = sum(org_volume),
            volume_perc = reduc_percentage(volume, org_volume)) %>%
  ungroup() %>%
  mutate(ietd = as.character(ietd))

data_plot2 <- BR_peak_volumes %>%
  filter(ietd == 2) %>%
  group_by(ietd, gi, year) %>%
  summarise(volume = sum(volume),
            org_volume = sum(org_volume),
            volume_perc = reduc_percentage(volume, org_volume)) %>%
  ungroup() %>%
  mutate(ietd = as.character(ietd)) %>%
  group_by(ietd, gi) %>%
  summarise(q10 = quantile(volume_perc, 0.1),
            q90 = quantile(volume_perc, 0.9))

data_plot <- data_plot %>%
  left_join(data_plot2, by = c("ietd", "gi"))

data_annual_rain <- data_frame(
  year = c(1:10),
  rain = (sapply(rains_list,sum))
)

data_plot3 <- BR_peak_volumes %>% # dot plots
  filter(ietd == 2) %>%
  group_by(ietd, gi, year) %>%
  summarise(volume = sum(volume),
            org_volume = sum(org_volume),
            volume_perc = reduc_percentage(volume, org_volume)) %>%
  ungroup() %>%
  filter(gi %% 0.5 == 0) %>%
  left_join(data_annual_rain, by = "year")

ggplot(data_plot) +
  geom_point(data = data_plot3, aes(gi, volume_perc, size = rain), alpha = 0.2, colour = "blue") +
  geom_ribbon(aes(x = gi, ymax = q10, ymin = q90), alpha = 0.4,show.legend = F) +
  scale_size(range = c(2,5)) +
  geom_line(data = data_plot, aes(gi, volume_perc)) +
  labs(y = "Target runoff volume reduction percentage [%]",
       x = "Required bioretention cell area [%]",
       size = "Rainfall depth [mm]") +
  coord_flip() +
  theme_bw(base_size = 7) +
  theme(legend.justification=c(0,1), legend.position=c(0,1),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))

ggsave("figures/figure7.pdf",  width = 90, height = 70, units = "mm")

save(BR_peak_volumes,data_plot,data_plot2, data_plot3, file='figures/Figure7.Rda')


library(tidyverse)
load("Figure7.Rda")

write.csv(data_plot, "fig7.csv")










