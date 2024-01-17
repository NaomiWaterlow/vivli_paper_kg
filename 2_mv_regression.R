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

####### Explore output for those with high index values: note only do MV on those with > 10000 samples so not all in there
store <- read_csv("mv_output/store_all.csv")
store_n <- store %>% mutate(combo = paste0(antibiotic,"_",bacteria))

# Bug-drug combinations with high index values
high_combos = as.data.frame(
  rbind(c("age_sex","levofloxacin_Staphylococcus aureus"),
        c("key_source","levofloxacin_Escherichia coli"),
        c("age_sex","levofloxacin_Klebsiella pneumoniae"),
        c("age_sex","levofloxacin_Pseudomonas aeruginosa"),
        c("age_sex","imipenem_Staphylococcus aureus"),
        c("age_sex","doripenem_Pseudomonas aeruginosa"),
        c("age_sex","cefepime_Escherichia coli"),
        c("age_sex","cefepime_Klebsiella pneumoniae")))

# Other high index possibilities
#c("age_sex","ceftolozane tazobactam_Klebsiella pneumoniae"),
#c("key_source","moxifloxacin_Staphylococcus aureus"),
#c("key_source","trimethoprim sulfa_Staphylococcus aureus"),
#c("key_source","levofloxacin_Staphylococcus aureus"),
#c("key_source","ciprofloxacin_Escherichia coli"),
#c("key_source","ceftaroline_Klebsiella pneumoniae"),
#c("key_source","ertapenem_Klebsiella pneumoniae"),
#c("key_source","ceftaroline_Pseudomonas aeruginosa")))

colnames(high_combos) <- c("grouping","combo")

store_nhc <- store_n %>% filter(combo %in% high_combos$combo)

length(unique(store_nhc$combo)) # 12
dim(high_combos) # 15
### Some high index not in MV as fewer than 10,000 samples
write.csv(store_nhc, "mv_output/store_nhighcombos.csv")
store_nhc %>% pivot_wider()

