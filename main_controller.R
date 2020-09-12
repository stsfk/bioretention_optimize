# Purpose:
#
# The script runs SWMM simulations of bioretentions of various surface areas under various rainfalls.
# 

# Libraries ---------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, parallel, stringr, swmmr)


# Set up clusters for SWMM simulation -----------------------------------

no_cores <- min(5, detectCores() - 1) # no_cores is the number of clusters for parallel simulation of SWMM
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
  
  file_location <- str_c(to_folder, "/",folder_prefix, folder_suffix,collapse = "")
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

write_br_input <- function(inp, br_area,br_width,drain_percent_br,width_imp,input_template, year, n){
  
  
  # change path of the rainfall data
  inp$timeseries$Value <- paste0("rain_", year, ".csv")
  
  # change bioretention implementation scenarios
  inp$lid_usage$Area[3] <- br_area
  inp$lid_usage$Width[3] <- br_width
  inp$lid_usage$FromImp[3] <- drain_percent_br
  
  # change the width of impervious area
  inp$subcatchments$Width <- width_imp
  
  # write input file
  write_inp(inp, paste0("./input_files/test",n,".inp"))
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
  dat <- read.table(input_file, skip = skip)
  dat <- matrix(unlist(dat[, 8]), ncol = 2, byrow = T) * conveter # convert to m3/min
  
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

grid <- expand.grid(year = 1:years, case = 1:cases) %>%
  tibble() %>%
  mutate(br_area = (area_imp * 0.02 + (case - 1) * (1.236 * 0.001))*43560, # surface area of bioretention cells, in ft2
         br_width = br_area/20, # width of of bioretention cells, assuming drainage path length = 20
         area_imp_lid = area_imp * 43560 - br_area, # area of impervious area excluding LID
         width_imp = area_imp_lid / 100, # width of imperious area excluding LID, assuming drainage path length = 100 
         drain_percent_br = min(br_area / area_imp_lid * br_AR * 100, 100)) # percent of impervious area drains to bioretention

for (i in 1:years){
  # evaluate bioretention implementation scenarios in different years
  
  # write input file
  dat <- grid %>%
    filter(year == i) %>%
    select(br_area,br_width,drain_percent_br,width_imp) %>% # change this line to the variable of interest
    map(as.character) %>%
    c(n = list(c(1:cases))) %>%
    pmap(write_br_input, inp = inp, year = i)
  
  # copy rainfall files to cluster
  for (k in 1:no_cores){
    rain_folder <- paste0(org_wd,"/data")
    move_file(paste0("rain_",i,".csv"), folder_suffix = NULL, file_suffix = i,
              file_prefix = "rain_", from_folder = rain_folder, to_folder = paths[[k]])
  }
  
  # create function to keep SWMM and rainfall files after execute a simulation
  remain_files <- c("swmm5.exe",paste0("rain",i,".csv"), "swmm5.dll")
  clean <- clean_after_run(remain_files)
  
  # path of the SWMM input file to be evaluated
  input_files <- list.files(path = "./input_files", pattern='^test.*\\.inp', recursive=TRUE)
  input_files <- paste("./input_files/test", seq_along(input_files),".inp",sep="")
  
  # start simulation
  le <- ceiling(length(input_files)/no_cores)
  
  for (j in 1:le){
    # move SWMM input file to clusters
    ind <- (j - 1)*no_cores + 1:no_cores
    file_2_copy <- as.list(paste0("test",ind,".inp")) 
    folder_suffix <- as.list(1:no_cores)
    file_suffix <- as.list(ind)
    arguments <- list(file_2_copy, folder_suffix, file_suffix) 
    arguments %>%
      pmap(move_file,"cluster","test", from_folder = "./input_files/")
    
    # run on clusters "cl" using run_swmm function
    parLapply(cl, paths, run_swmm) 
    
    # move file
    for (k in 1:no_cores){
      move_file("outflow.txt",folder_suffix = i, file_suffix = ind[k], folder_prefix = "year",
                file_prefix = "br",from_folder = paths[[k]], to_folder = "./results")
      move_file("test1.rpt",folder_suffix = i, file_suffix = ind[k], folder_prefix = "year",
                file_prefix = "br",from_folder = paths[[k]], to_folder = "./results")
    }
    
    # clean cluster for evaluating new input files
    lapply(paths,clean)
  }
  
  # clean cluster
  
  delete_file(input_files)
  clean <- clean_after_run(c("swmm5.exe"))
  lapply(paths,clean)
  
  # write surf and ud series files
  file_location <- str_c(org_wd, "/", "year", i, "/", collapse = "")
  setwd(file_location)
  
  br_files <- list.files(pattern = '^br.*\\.txt')
  br_files <- paste0("br", c(1:length(br_files)), ".txt")
  
  for (j in seq_along(br_files)) {
    input_file = br_files[j]
    output_file1 = paste0("surf",j,".txt")
    output_file2 = paste0("ud",j,".txt")
    write_series(input_file, output_file1, output_file2)
  }
  
  delete_file(br_files)  
  setwd(org_wd)
}

stopCluster(cl)

