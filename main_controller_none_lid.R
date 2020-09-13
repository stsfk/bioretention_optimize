# Purpose:
#
# The script runs SWMM simulations of the impervious catchment under various rainfalls.
# 

# Libraries ---------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, parallel, stringr, swmmr)


# Set up clusters for SWMM simulation -----------------------------------

no_cores <- min(10, detectCores() - 1) # no_cores is the number of clusters for parallel simulation of SWMM
paths <- as.list(str_c(getwd(),"/cluster",c(1:no_cores)))
cl <- makeCluster(no_cores) # start cluster
org_wd <- getwd()

# Functions ---------------------------------------------------------------

delete_file <- function(x) {
  # This function is used to delete file with name x
  if(!identical(x, character(0))){
    if (any(file.exists(x))) {
      file.remove(x)
      T
    } else {
      F
    }
  } else {
    F
  }
}

run_swmm <-  function(path, pattern = '^test.*\\.inp') {
  # This function is used to run SWMM inside path folder, this allows SWMM simulation to be parallelized over different clusters 
  # The name of the input file follows pattern, '^test.*\\.inp' runs SWMM input file started with "test"
  
  # Input:
  #   file_name is the name of SWMM input file
  #   path is folder address
  org_wd <- getwd()
  
  setwd(path)
  file_name = list.files(pattern = pattern)
  
  if (!is.null(file_name)) {
    file_name <- file_name[1]
    com <- paste("swmm5.exe", file_name, "test1.rpt test1.out", sep = " ")
    system(com)
    setwd(org_wd)
    return(T)
  } else {
    setwd(org_wd)
    return(f)
  }
}

move_file <- function(file_name, folder_suffix, file_suffix, 
                      folder_prefix = NULL, file_prefix = NULL, 
                      from_folder = getwd(), to_folder = getwd()){
  # This function is to move output to a give output location
  # Input:
  #   file_name: name of the file to be moved
  #   folder_suffix: suffix of folder, e.g. a number
  #   file_suffix: suffix of file, e.g. a number 
  #   file_preffix: prefix of a file, e.g. "BC"
  #   folder_sufix: folder_sufix of the folder, e.g. some category, "year"
  #   folder_prefix: prefix of file, e.g. "year"
  # Output:
  #   is_moved: T or F
  org_wd <- getwd()
  
  setwd(from_folder)
  
  file_location <- str_c(to_folder, "/",folder_prefix, folder_suffix, collapse = "")
  if (!dir.exists(file_location)) {
    dir.create(file_location)
  }
  
  pattern <- "\\."
  loc <- unlist(str_locate(file_name,pattern))[,2]
  to_name <- str_c(file_location, "/",
                   file_prefix,file_suffix,
                   str_sub(file_name, loc,-1),
                   collapse = T)
  
  file.copy(file_name,to_name,overwrite = T)  
  
  setwd(org_wd)
}

write_imp_input <- function(inp, year){
  
  
  # change path of the rainfall data
  inp$timeseries$Value <- paste0("rain_", year, ".csv")
  
  # write input file
  write_inp(inp, paste0("./input_files/test", year,".inp"))
}

clean_after_run <- function(x) {
  # This is a function factory
  # that produces functions to clean folder, only keep file "x"
  function(path){
    org_wd <- getwd()
    setwd(path) 
    delete_file(setdiff(list.files(path = path),x))
    setwd(org_wd)
  }
} 

write_series <- function(input_file, output_file1, output_file2, skip = 9, conveter = 1.699011){
  # this function is to write time series files of underdrain flow and surface runoff
  # the resolution is 1-minute
  dat <- read.table(input_file, skip = skip)
  dat <- matrix(unlist(dat[, 8]), ncol = 2, byrow = T) * conveter # convert from CFS to m3/min
  
  file_name_surf <- str_c(output_file1, collapse = "")
  file_name_ud <- str_c(output_file2, collapse = "")
  
  write.table(
    dat[, 1],
    file  = file_name_surf,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE ,
    sep = ""
  )
  write.table(
    dat[, 2],
    file  = file_name_ud,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE ,
    sep = "")
}


# constant ----------------------------------------------------------------
years <- 10 # 10 years of rainfalls
cases <- 81 # number of bioretention implementation scenarios, surface area ranges from 2% to 10% of the contributing drainage area
area_imp <- 1.236 # area of catchment, 5,000 m2

inp <- read_inp("./data/catchment_bioretention.inp") # read in template

br_AR = 50 # bioretention cells' maximum contributing drainage area capture ratio, i.e., it treats area of up to 50 times of its surface area

# main --------------------------------------------------------------------

# prepare variables to be changed for evaluating bioretention performance
years <- 1:10

# start SWMM evaluation

# write input file for each year
sapply(1:10, write_imp_input, inp = inp)

# copy rainfall files to cluster
file2copy <- paste0(org_wd, "/data/rain_", 1:10, ".csv")
file2 <- paste0(org_wd, "/cluster", 1:10, "/rain_", 1:10, ".csv")
file.copy(file2copy, file2, overwrite = T)

# copy swmm input files to cluster
file2copy <- paste0(org_wd, "/input_files/test", 1:10, ".inp")
file2 <- paste0(org_wd, "/cluster", 1:10, "/test", 1:10, ".inp")
file.copy(file2copy, file2, overwrite = T)

# run swmm simulations on clusters "cl" using run_swmm function
parLapply(cl, paths, run_swmm)

# copy swmm results

file2copy <- paste0(org_wd, "/cluster", 1:10, "/outflow.txt")
file2 <- paste0(org_wd, "/results/year", 1:10, "/br0.txt")
file.copy(file2copy, file2, overwrite = T)

file2copy <- paste0(org_wd, "/cluster", 1:10, "/test1.rpt")
file2 <- paste0(org_wd, "/results/year", 1:10, "/br0.rpt")
file.copy(file2copy, file2, overwrite = T)

# clean clusters
remain_files <- c("swmm5.exe", "swmm5.dll")
clean <- clean_after_run(remain_files)

lapply(paths, clean)

# write surface (surf)) and uderdrain (ud) series files for each year
for (i in 1:10){
  file_location <- str_c(org_wd, "/results/", "year", i, "/", collapse = "")
  setwd(file_location)
  
  output_file1 = paste0("surf0.txt")
  output_file2 = paste0("ud0.txt")
  write_series("br0.txt", output_file1, output_file2)
  
  delete_file("br0.txt")  
  setwd(org_wd)
}

# stop clusters
stopCluster(cl)
