library(tidyverse)
library(modelr)
library(stringr)

# Read data ---------------------------------------------------------------
data <- read.table("clipboard", header = T)

data <- data %>%
  select(gi, rv, ff, q10_2, q10_4, w_2, w_4
  )

# Functions ---------------------------------------------------------------
predict_fac <- function(model){
  # function factory for prediction functions
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

q10_2_m <- loess(q10_2 ~ gi, data = data)
q10_2_f <- predict_fac(q10_2_m)

q10_4_m <- loess(q10_4 ~ gi, data = data)
q10_4_f <- predict_fac(q10_4_m)

w_2_m <- loess(w_2 ~ gi, data = data)
w_2_f <- predict_fac(w_2_m)

w_4_m <- loess(w_4 ~ gi, data = data)
w_4_f <- predict_fac(w_4_m)

f <- function(FF, PF, RV){
  # A function factory for calculating system-wide score
  function(D, A, a, b, case){
    # case is reduction coefficient case
    ratio <- D/100
    if (case == 1){
      coe <- (0.5 + ratio*5/6)
    } else if (case == 2){
      coe <- 1
    } else {
      coe <- (1.5 - ratio*5/6)
    }
    
    Ap <- A/D*100 # area percentage
    p_FF <- cummin(FF(Ap)) # performance increase with increased area of BC
    p_PF <- cummin(PF(Ap))
    p_RV <- cummin(RV(Ap))
    
    temp <- D*(p_FF + a*p_PF*coe + b*p_RV)
    
    return(temp)
  }
}

find_max_f <- function(f_function, a, b){
  x <- 2 + 0.1*(0:80) # avaiable BC area
  control_area1 <- rep(0, length(x)) # optimal control area
  control_area2 <- rep(0, length(x)) 
  control_area3 <- rep(0, length(x)) 
  
  for (i in seq_along(x)){
    D <- seq(from = x[i]*10, 100, by = 0.1) # possiable solutions
    
    score1 <- f_function(D,  x[i], a, b, 1) # case 1
    score2 <- f_function(D,  x[i], a, b, 2) # case 2
    score3 <- f_function(D,  x[i], a, b, 3) # case 3
    
    control_area1[i] <- D[which.max(score1)]
    control_area2[i] <- D[which.max(score2)]
    control_area3[i] <- D[which.max(score3)]
  }
  
  return(cbind(control_area1, control_area2, control_area3))
}

f_2 <- f(ff_f, w_2_f, rv_f)
f_4 <- f(ff_f, q10_4_f, rv_f)

para <- expand.grid(c(0,0.5,1,2,5),c(0,0.5,1,2,5))

simu_fun <- function(f_x, a, b){
  q <- find_max_f(f_x, a, b)
  q <- as.data.frame(q)
  q$a <- a
  q$b <- b
  q$gi <- 2 + 0.1*(0:80)
  q <- q %>%
    gather(type, value, -a, -b, -gi)
  
  return(q)
}

results <- vector("list", length(para[[1]]))
for (i in 1:length(para[[1]])){
  para_a <- as.numeric(para[i,1])
  para_b <- as.numeric(para[i,2])
  
  result_2 <- simu_fun(f_2,para_a,para_b)
  result_2$storm <- "[25,50)"
  result_4 <- simu_fun(f_4,para_a,para_b)
  result_4$storm <- "[75,100)"
  
  results[[i]] <- rbind(result_2,result_4)
}

results <- bind_rows(results)

# plot --------------------------------------------------------------------
results$a <- as.factor(results$a)
levels(results$a) <- c("a = 0",
                       "a = 0.5",
                       "a = 1",
                       "a = 2",
                       "a = 5")

results$b <- as.factor(results$b)
levels(results$b) <- c("b = 0",
                       "b = 0.5",
                       "b = 1",
                       "b = 2",
                       "b = 5")

results$type <- as.factor(results$type)
levels(results$type) <- c("pro-extensive", "neutral", "pro-intensive")

results <- results %>% 
  filter(a != "a = 0.5",
         b != "b = 0.5")

ggplot(results, aes(gi, value, colour = type, linetype = storm)) +
  geom_line(size = 0.5) +
  facet_grid(b~a) +
  labs(x = "Available bioretention cell area as the percentage of catchment area [%]",
       y = "Optimal CDA as the percentage of the catchment area [%]",
       linetype = "Storm groups",
       colour = "Preference") +
  theme_bw(base_size = 7)

ggsave("Figure9.pdf", width = 190, height = 120, units = "mm")

save(results, file='figure9.Rda')



