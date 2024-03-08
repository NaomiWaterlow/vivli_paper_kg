##### Run analysis

##### MICAG Tool for screening and plotting MIC by sub_group
library(data.table);library(ggplot2);library(cowplot); library(dplyr)
theme_set(theme_bw())
### Load in the functions
source("1_functions.R")

# read in the data
# option to load in own data here. Must be same format. 
full_data <- as.data.table(read.csv("data/full_data.csv"))

# specify which bugs are of interest
bacteria_to_use <- unique(full_data$organism_clean) 

## What characteristic to look at. (Note: Must match column name)
characteristics <- c("age_group", "key_source", "country", "income_grp", "who_region") # c("age_group", "key_source") #Can run additional options: 
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

## for who region
index_who_results <- read.csv("plots/gender_who_regionindex_store.csv")

summarised_n_index_who <- index_who_results %>%
  group_by(antibiotic, organism_clean, gender) %>%
  summarize(Sum_N = sum(N))

temp_n_max_index_who <- index_who_results %>%
  group_by(antibiotic, organism_clean, gender) %>%
  summarize(max_index = max(dff))

summarised_n_index_who <- cbind(summarised_n_index_who, temp_n_max_index_who[, -c(1:3)])
rm(temp_n_max_index_who)

summarised_n_index_who$grouping <- "who"


n_index_df <- rbind(summarised_n_index_age, summarised_n_index_source, summarised_n_index_who)

rm(summarised_n_index_age, summarised_n_index_source, summarised_n_index_who)

## Calculating correlation coefficient
correlation_result <- n_index_df %>%
  group_by(grouping, organism_clean, gender) %>%
  summarize(correlation_coefficient = cor(Sum_N, max_index, use = "complete.obs"))

correlation_result$vjust <- if_else(correlation_result$gender == "f", -1, 1)

## Labels for plot
supp.lab <- c("Age", "Source", "WHO region")
names(supp.lab) <- c("age", "source","who")

## Plot to explore if correlation
## Remove low numbers
table(n_index_df$Sum_N)
colnames(n_index_df) <- c("antibiotic", "organism_clean","Gender","Sum_N","max_index","grouping")
colnames(correlation_result) <- c("grouping", "organism_clean","Gender","correlation_coefficient","vjust")

ggplot(n_index_df %>% filter(Sum_N > 100), aes(x= Sum_N, y = max_index, color = Gender, shape=Gender)) +
  facet_grid(grouping ~ organism_clean, labeller = labeller(grouping = supp.lab)) +
  geom_point() +
  geom_smooth(method='lm',aes(fill=Gender))+
  scale_color_manual(values=c("#7FB069","#805D93")) +
  scale_fill_manual(values=c("#7FB069","#805D93")) +
  scale_shape_manual(values=c(16,1)) +
  scale_x_continuous(limits = c(0,100000))+
  geom_text(data = correlation_result,
            aes(x = 90000, y = 0.69, label = paste("R = ", round(correlation_coefficient, digits = 3)), vjust = vjust),
            size = 2
            ,show.legend = FALSE
            ) +
  xlab("Number of samples") +
  ylab("Maximum difference in MIC across groupings") +
  theme(strip.text.x = element_text(face = "italic"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))

ggsave(paste0("plots/", "figure2_indexgroupings.pdf"), height = 7, width = 10)
  