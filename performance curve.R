# Figure 8 ----------------------------------------------------------------

library(tidyverse)
library(modelr)
library(stringr)

# Start -------------------------------------------------------------------

data <- read.table("clipboard", header = T)

data[data > 100] <- 100 # deal with the prediction errors

predict_fac <- function(model){
  # function factory for generate predictions
  function(x){
    df <- data.frame(gi = x)
    out <- unlist(add_predictions(df, model)[,2])
    out[out > 100] <- 100
    
    out[out < 0] <- 0
    
    return(out)
  }
}

rv_m <- loess(rv ~ gi, data = data)
rv_f <- predict_fac(rv_m)
ff_m <- loess(ff ~ gi, data = data)
ff_f <- predict_fac(ff_m)

find_requr_bc <- function(tar_pec, fun){
  # find required BC area
  grids <- seq(from = 2, to = 10, by = 0.0001)
  eva <- fun(grids)
  
  ind <- which.min(abs(tar_pec - eva))
  grids[ind]
}

FF_area_req <- sapply(c(7:9)*10, find_requr_bc, fun = ff_f)
FF_df <- data_frame(FF_thd = c(7:9)*10,
                    FF_area_req = FF_area_req)
names(FF_area_req) <- paste0("FF ",(c(7:9)*10), "%")

RV_area_req <- sapply(c(5:7)*10, find_requr_bc, fun = rv_f)
RV_df <- data_frame(RV_thd = c(5:7)*10,
                    RV_area_req = RV_area_req)
names(RV_area_req) <- paste0("RV ",(c(5:7)*10), "%")

req <- c(FF_area_req, RV_area_req)
req <- as.data.frame(matrix(unname(req), ncol = 2, byrow = F))
names(req) <- c("FF", "RV")

# process -----------------------------------------------------------------
loess_model <- function(data) {
  loess(value ~ gi, data = data)
}

data_process <- data %>%
  select(-ff, -rv) %>%
  gather(item, value, -gi) %>%
  group_by(item) %>%
  nest() %>%
  mutate(model = map(data, loess_model),
         pred = map2(data, model, add_predictions)
  ) %>%
  unnest(pred)

data <- data_process %>%
  select(-value) %>%
  spread(item, pred)

temp <- data[1,] # this is to plot horizontal lines on the left
temp[1,-1] <- 0
data <- rbind(temp, data)
data[data > 100] <- 100
data[data < 0] <- 0

data_process <- data %>%
  gather(item, value, -gi) %>%
  group_by(item)%>%
  mutate(value = cummax(value)) %>%
  ungroup()

data_process <- data_process %>%
  group_by(item) %>%
  mutate(id = 1:n(),
         max_ind = which.max(value),
         keep = id <= max_ind) %>%
  filter(keep == T) %>%
  select(gi, item, value) %>%
  ungroup() 

ff_names <- paste0("FF reduction ",(c(7:9)*10), "%")
rv_names <- paste0("RV reduction ",(c(5:7)*10), "%")

result <- vector("list", 9)
k <- 1
for (i in 1:3){
  for (j in 1:3){
    rq <- max(req[i,1],req[j,2])
    
    df <- data_process
    temp <- df$gi
    temp[temp < rq] <- rq
    df$gi <- temp
    
    df <- df %>%
      filter(value < 100)
    
    df$ff <- ff_names[[i]]
    df$rv <- rv_names[[j]]
    
    result[[k]] <- df
    
    k <- k + 1
  }
}

result <- bind_rows(result)

result <- result %>%
  mutate(perc = str_sub(item, -1, -1),
         type = str_sub(item, 1, 1))

result$perc <- as.factor(result$perc)
levels(result$perc) <- c("[0,25)", "[25,50)","[50,75)","[75,100]")

result$type <- as.factor(result$type)
levels(result$type) <- c("q10","q.weighted")

ggplot(result, aes(value, gi, linetype = perc,  colour = type)) +
  geom_line(size = 0.4) +
  scale_linetype_manual(values = c("solid","dotted","dotdash","dashed")) +
  scale_y_continuous(limits = c(3.5,10), breaks = c(2:5)*2) +
  facet_grid(ff~rv) +
  scale_x_continuous(limits = c(0,100)) +
  theme_bw(base_size = 7) +
  theme(legend.key.width = unit(0.7, "cm")) +
  labs(x = "Target peak flow reduction percentage [%]",
       y = "Required bioretention cell area [%]",
       linetype = "Storm groups",
       colour = "Curve types")

ggsave("figure8.pdf", width = 190, height = 120, units = "mm")

save(result, file='Figure8.Rda')