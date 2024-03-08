##### Run analysis

##### MICAG Tool for screening and plotting MIC by sub_group
library(data.table);library(ggplot2);library(cowplot)

### Load in the functions
source("1_functions.R")

# read in the data
# option to load in own data here. Must be same format. 
full_data <- as.data.table(read.csv("data/full_data.csv"))

# specify which bugs are of interest
bacteria_to_use <- unique(full_data$organism_clean) 

## What characteristic to look at. (Note: Must match column name)
characteristics <- c("age_group", "key_source") #Can run additional options: 
#"key_source" # "age_group" # country # income_grp #who_region
#

### Run initial plot and index generation
plot_generation_MICAG(full_data, bacteria_to_use, characteristics, gender_options = c("F","T")) 

### Run by time 
plot_generation_bytime_MICAG(full_data, bacteria_to_use, characteristics, gender_options = c("F","T")) 
