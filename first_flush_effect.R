## FFE reduction calculation


# Libraries ---------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, xts, zoo, lubridate)

# Set up ------------------------------------------------------------------

org_wd <- getwd()

# Identify storm starts and ends ------------------------------------------

read_rain <- function(file){
  # This function read rainfall time series in "data" folder, 
  # the output is a vector of rainfall at different time steps
  
  temp <- readLines(file)
  
  rain_depths <- as.numeric(str_sub(temp, 18, -1))
  rain_depths<- rain_depths/6*25.4 # convert from inch/h to mm per 10 min
  
  return(c(rain_depths))
}

sta_end <- function(x, par.inter) {
  # This function is used to find start and end time of a storm event based on the dry spell threshold, 
  # i.e., inter-event time definition (IETD).
  # Input:
  #   x is time series
  #   par.inter is the dry spell interval
  # Output:
  #   sta_end_data is a list containing the start end information, i.e., the index of the rainfall event in a time series
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
    # in SWMM,  rainfall intensity is assumed to last for the coming time step, ajusted to intensity of last time step
    rains_list[[i]] <- c(0,read_rain(file))
  }
}

inter <- as.list((1:4)*48) # IETD = 8,16,24,32 hour
sta_ends <- list() # sta_ends stores the start and end index for different par.inter
for (i in seq_along(rains_list)){
  rain <- c(rains_list[[i]])
  sta_ends[[i]] <- lapply(inter, sta_end, x = rain)
}

# Read surface runoff time series -----------------------------------------
setwd(org_wd)

x_gis <- list()
x_gis_aggre <- list() # x_gis_aggre[[j]][[i]], j different years, i different BC %
x_org_aggre <- list() # orginal runoff without BC cover
for (j in 1:10){
  file_location <- str_c(org_wd, "/results/year",j,"/",collapse = "")
  files <- list.files(pattern='surf.*\\.txt', recursive=TRUE, path = file_location)
  files <- paste("surf", c(1:(length(files))),".txt",sep="")
  x_gis <- list()# x_gis is a list containing GI of different % area, 2:9 in this case
  
  for (i in 1:81) {
    fid <- i
    file.name <- paste(file_location, files[fid], sep = "")
    x_gis[[i]] <- unname(unlist(read.table(file.name, header = F)))
  }
  org_name <- paste(file_location,"surf0.txt",sep = "")
  x_org_aggre[[j]] <- unname(unlist(read.table(org_name)))
  x_gis_aggre[[j]] <- x_gis
}

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

# FF calculation ----------------------------------------------------------
depths <- c(1:6)*25 # 0.5 cm to 3 cm

first_true <- function (x) {
  # Return index of first TRUE
  # If all FALSE, return length of the the vector
  if (sum(x) >= 1) {
    return (min(which(x == TRUE)))
  } else {
    return(length(x))
  }
}

find_FF_id_depth <- function(x, breaks, depth){
  # This function is to find the id and depth of the initial portion
  labels <- 1:(length(breaks) - 1)
  
  x <- data.frame(x = x, id = seq_along(x)) %>%
    mutate(interval = cut(id, breaks = breaks, include.lowest = TRUE, right = F, labels = labels),
           interval = as.numeric(interval)) %>%
    group_by(interval) %>%
    mutate(cum_x = cumsum(x),
           FF_id = (cum_x >= depth)) %>%
    summarise(FF_id = first_true(FF_id),
              FF_depth = cum_x[FF_id],
              ratio = ifelse(FF_depth > depth, (FF_depth - depth)/(x[FF_id]), 0)) # for cases last runoff larger than depth thd
  return(x)
}

find_FF_id_depths <- function(x_org_aggre, x_intervals_aggre, depth){
  org_runoff_depths <- vector("list", 10)
  for (i in 1:10){
    x <- x_org_aggre[[i]]
    org_runoff_depth <- vector("list", length(inter)) # interval
    for (j in 1:length(inter)) {
      breaks <- x_intervals_aggre[[i]][[j]]
      temp <- find_FF_id_depth(x, breaks, depth)
      temp$interval <- as.numeric(temp$interval)
      temp$year <- i
      temp$ietd <- j
      org_runoff_depth[[j]] <- temp
    }
    org_runoff_depth <- bind_rows(org_runoff_depth)
    org_runoff_depths[[i]] <- org_runoff_depth
  }
  temp <- bind_rows(org_runoff_depths)
  temp$depth <- depth
  return(temp)
}

temp <- lapply(depths, find_FF_id_depths, x_org_aggre = x_org_aggre, x_intervals_aggre = x_intervals_aggre)
org_depths <- bind_rows(temp)

find_FF_depth_gi <- function(x_aggre, org_depths, year_ind, ietd_ind, depth_ind){
  x <- x_aggre[[year_ind]]
  breaks <- x_intervals_aggre[[year_ind]][[ietd_ind]]
  labels <- 1:(length(breaks) - 1)
  org_depth <- org_depths %>%
    dplyr::filter(year == year_ind, ietd == ietd_ind, depth == depth_ind)
  
  x <- data.frame(x = x, id = seq_along(x)) %>%
    mutate(interval = cut(id, breaks=breaks, include.lowest = TRUE, right = F, labels = labels),
           interval = as.numeric(interval)) %>%
    left_join(org_depth, by = "interval") %>%
    select(id, x, interval, FF_id, FF_depth, ratio) %>%
    group_by(interval) %>%
    mutate(cum_x = cumsum(x)) %>%
    summarise(FF_depth_GI = cum_x[FF_id [1]] - x[FF_id [1]]*ratio[1],
              FF_depth_org = min(mean(FF_depth), depth_ind))
  return(colSums(x)[2:3])
}

FF_depths <- vector("list", 81)
for (i in 1:81){
  # loop for each gi area
  x_aggre <- vector("list", 10)
  for (j in 1:10){
    x_aggre[[j]] <- x_gis_aggre[[j]][[i]] # i is gi area, j is year
  }
  
  df1 <- expand.grid(year_ind = c(1:10), ietd_ind = c(1:length(inter)), depth_ind = depths)
  df2 <- as.list(df1) %>%
    pmap(find_FF_depth_gi, x_aggre = x_aggre, org_depths = org_depths)
  df2 <- data.frame(matrix(unlist(df2), ncol = 2, byrow = T))
  
  df <- cbind(df1, df2)
  df <- df %>%
    rename(gi = X1, org = X2) %>%
    mutate(perc = gi/org*100) %>%
    group_by(ietd_ind, depth_ind) %>%
    summarise(sum_gi = sum(gi)*2/10, # convert to per ha per year
              sum_org = sum(org)*2/10,
              mean = sum(gi)/sum(org)*100,
              q10 = quantile(perc, 0.1),
              q90 = quantile(perc, 0.9))
  
  df$gi <- i * 0.1 + 1.9
  
  FF_depths[[i]] <- df
}

FF_depths <- bind_rows(FF_depths)

# Figure 4 ----------------------------------------------------------------

data_plot <- FF_depths

data_plot$depth_ind <- as.factor(data_plot$depth_ind)
levels(data_plot$depth_ind) <- c("0.5 cm","1.0 cm","1.5 cm","2.0 cm","2.5 cm", "3.0 cm")

data_plot$ietd_ind <- as.factor(data_plot$ietd_ind)
levels(data_plot$ietd_ind) <- c("IETD = 8h","IETD = 16h","IETD = 24h","IETD = 32h")

data_analysis <- data_plot %>%
  filter(depth_ind == "2.5 cm", ietd_ind == "IETD = 24h") %>%
  select(gi, sum_gi, sum_org) %>%
  View()

data_plot <- FF_depths

data_plot2 <- data_plot %>%
  mutate(depth_ind = as.numeric(as.character(depth_ind)),
         tar_percentage = 100 - mean,
         tar_q10 = 100 - q10,
         tar_q90 = 100 - q90) %>%
  filter(depth_ind %% 25 == 0,
         depth_ind < 175)

data_plot2$depth_ind <- as.factor(data_plot2$depth_ind)
levels(data_plot2$depth_ind) <- c("0.5 cm","1.0 cm","1.5 cm","2.0 cm","2.5 cm", "3.0 cm")

data_plot2$ietd_ind <- as.factor(data_plot2$ietd_ind)
levels(data_plot2$ietd_ind) <- c("IETD = 8h","IETD = 16h","IETD = 24h","IETD = 32h")


ggplot(data_plot2, aes(gi, tar_percentage)) +
  geom_ribbon(aes(x = gi, ymax = tar_q10, ymin = tar_q90,fill = depth_ind),alpha = 0.2,show.legend = F) +
  geom_line(aes(colour = depth_ind,linetype = depth_ind), size = 0.5) +
  facet_wrap(~ietd_ind) +
  scale_y_continuous(limits = c(40,100), breaks = c(0:3)*20+40) +
  labs(y = "Target reduction percentage of first flush volume [%]",
       x = "Required bioretention cell area [%]",
       colour = "Initial runoff",
       linetype = "Initial runoff") +
  coord_flip() +
  theme_bw() +
  theme(legend.key.width = unit(1, "cm")) 

ggsave("./figures/Figure4.pdf", width = 190, height = 120, units = "mm") 

save(FF_depths, file='./figures/Figure4.Rda')

# Analysis ----------------------------------------------------------------

data_plot <- FF_depths

data_plot2 <- data_plot %>%
  mutate(depth_ind = as.numeric(as.character(depth_ind)),
         tar_percentage = 100 - mean,
         tar_q10 = 100 - q10,
         tar_q90 = 100 - q90) %>%
  filter(depth_ind %% 25 == 0,
         depth_ind < 175)

data_analysis <- data_plot2 %>%
  filter(ietd_ind == 2,
         depth_ind == 125) %>%
  select(gi,tar_percentage)

write.csv(data_analysis, "./figures/fig4.csv")


