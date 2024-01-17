##### MICAG regression analysis for specified drugs
# Specify drug / bug combination to run regression analysis
library(data.table); library(MASS); library(tidyverse)

# read in data (can use alternative)
full_data_orig <- as.data.table(read.csv("data/full_data.csv"))

# Create folder for outputs
if(!file.exists("mv_output")){dir.create("mv_output")} 

# Load function
source("1_functions.R")

######*********************** SPECIFY WHICH BUG + DRUG ************************#################
# specify which bacteria and antibiotic of interest
# NOTE - must match spelling in the data table
mv_analysis(full_data_orig, "levofloxacin", "Staphylococcus aureus")

######*********************** OR CYCLE THROUGH ALL ************************#################
top_bacteria <- unique(full_data_orig$organism_clean)
store <- c() # where will store outputs

### Run through all 
for(i in top_bacteria){
  target_bug <- i
  full_data_bug <- full_data_orig %>% filter(organism_clean == target_bug)
  # which resistances were tested for? 
  t <- full_data_bug %>% group_by(antibiotic) %>% summarise(n=n()) 
  # Want more than 10000 samples or > 10%
  top_drugs <- t %>% filter(n > 10000)
  for(j in top_drugs$antibiotic){
    target_antibiotic <- j
    m <- mv_analysis(full_data_orig, target_antibiotic, target_bug)
    print(paste0("Run for ",target_antibiotic, " ", target_bug))
    ## Store
    store <- rbind(store, m %>% mutate(antibiotic = target_antibiotic, 
                                       bacteria = target_bug))
  }
}

write.csv(store, "mv_output/store_all.csv")
  
