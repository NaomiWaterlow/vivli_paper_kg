##### Run analysis

##### MICAG Tool for screening and plotting MIC by sub_group
library(data.table);library(ggplot2);library(cowplot); library(dplyr)

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
plot_generation_MICAG(full_data, bacteria_to_use, characteristics) 

### Run by time 
plot_generation_bytime_MICAG(full_data, bacteria_to_use, characteristics) 




## We would like to see the effect of sample size on max index value, to see whether over
## or undersampling affects whether a large max index value is observed.

## This is therefore code to generate regression plots of n vs max index for each grouping variable,
## separated by bug.

index_age_results <- read.csv("plots/gender_age_groupindex_store.csv")


summarised_n_index_age <- index_age_results %>%
  group_by(antibiotic, organism_clean, gender) %>%
  summarize(Sum_N = sum(N))

temp_n_max_index_age <- index_age_results %>%
  group_by(antibiotic, organism_clean, gender) %>%
  summarize(max_index = max(dff))

summarised_n_index_age <- cbind(summarised_n_index_age, temp_n_max_index_age[, -c(1:3)])
rm(temp_n_max_index_age)

summarised_n_index_age$grouping <- "age"

## for source
index_source_results <- read.csv("plots/gender_key_sourceindex_store.csv")


summarised_n_index_source <- index_source_results %>%
  group_by(antibiotic, organism_clean, gender) %>%
  summarize(Sum_N = sum(N))

temp_n_max_index_source <- index_source_results %>%
  group_by(antibiotic, organism_clean, gender) %>%
  summarize(max_index = max(dff))

summarised_n_index_source <- cbind(summarised_n_index_source, temp_n_max_index_source[, -c(1:3)])
rm(temp_n_max_index_source)

summarised_n_index_source$grouping <- "source"


n_index_df <- rbind(summarised_n_index_age, summarised_n_index_source)

rm(summarised_n_index_age, summarised_n_index_source)

## Plot

ggplot(n_index_df, aes(x= max_index, y = Sum_N, color = gender))+
  facet_wrap(~ grouping + organism_clean, ncol = 4)+
  geom_point()
  
