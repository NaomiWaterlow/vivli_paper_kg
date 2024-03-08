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
store_nhc %>% pivot_wider(id_cols = "parameter")


## The code below was designed to generate figures for the nhighcombos output
# Read in data of interest
mv_df <- read.csv("mv_output/store_nhighcombos.csv")

# Split by antibiotic first, then generate dataframes for each subgrouping category, as figures will first be 
# split by antibiotic, then by grouping category, then each bacterial species will be plotted on the same
# graphs using different colours.

## Grab the first 3 letters of the parameters column values and put into a new column named "grouping"
mv_df$grouping <- substr(mv_df$parameter, 1, 3)

mv_df$sig <- ifelse(mv_df$p.value < 0.001, "***", ifelse(mv_df$p.value<0.01, "**", ifelse(mv_df$p.value<0.05, "*", "")))

## Split the mv_df into the different antibiotics
antibiotic_split <- split(mv_df, mv_df$antibiotic)

## Split each of the antibiotic dataframes by the grouping value
final_dataframes <- lapply(antibiotic_split, function(sub_df) {
  split(sub_df, sub_df$grouping)
})

# Convert each dataframe in the list to separate objects in the global environment
for (antibiotic_name in names(final_dataframes)) {
  for (group_name in names(final_dataframes[[antibiotic_name]])) {
    obj_name <- paste(antibiotic_name, group_name, sep = "_")
    assign(obj_name, final_dataframes[[antibiotic_name]][[group_name]], envir = .GlobalEnv)
  }
}


## Remove the grouping variable from each parameter value:
levofloxacin_age$parameter <- sub("^age_group", "", levofloxacin_age$parameter)
cefepime_age$parameter <- sub("^age_group", "", cefepime_age$parameter)
doripenem_age$parameter <- sub("^age_group", "", doripenem_age$parameter)
levofloxacin_gen$parameter <- sub("^gender", "", levofloxacin_gen$parameter)
cefepime_gen$parameter <- sub("^gender", "", cefepime_gen$parameter)
doripenem_gen$parameter <- sub("^gender", "", doripenem_gen$parameter)
levofloxacin_key$parameter <- sub("^key_source", "", levofloxacin_key$parameter)
cefepime_key$parameter <- sub("^key_source", "", cefepime_key$parameter)
doripenem_key$parameter <- sub("^key_source", "", doripenem_key$parameter)
levofloxacin_who$parameter <- sub("^who_region", "", levofloxacin_who$parameter)
cefepime_who$parameter <- sub("^who_region", "", cefepime_who$parameter)
doripenem_who$parameter <- sub("^who_region", "", doripenem_who$parameter)

require("ggplot2")
require(scales)
install.packages("ggstance")
require("ggstance")
require(viridis)

levo_age <- ggplot(levofloxacin_age, aes(x= parameter, y=Odds, color = bacteria, group = bacteria))+
  geom_hline(yintercept = 1, color = "black")+
  geom_line(size = 2)+
  geom_point(size = 4)+
  theme_bw()+
  scale_y_continuous(expand= c(0,0), limits = c(0,2.5))+
  xlab("Age Group")+
  labs(color = "Species")+
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.text = element_text(size = 12))+
  scale_color_viridis_d(option = "D", labels = label_wrap(10))+
  geom_text(
    data = levofloxacin_age,
    aes(x= parameter, y=max(Odds)+0.2, color = bacteria, label = sig),
    position = position_dodgev(height = 0.2), size = 7, show.legend = FALSE)+
  guides(color = guide_legend(keyheight = 2.5))

levo_gen <- ggplot(levofloxacin_gen, aes(x= parameter, y=Odds, color = bacteria, group = bacteria))+
  geom_hline(yintercept = 1, color = "black")+
  geom_point(size = 4)+
  theme_bw()+
  scale_y_continuous(expand= c(0,0), limits = c(0,1.7))+
  xlab("Gender")+
  labs(color = "Species")+
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.position = "none")+
  scale_color_viridis_d(option = "D", labels = label_wrap(10))+
  geom_text(
    data = levofloxacin_gen,
    aes(x= parameter, y=max(Odds)+0.2, color = bacteria, label = sig),
    position = position_dodgev(height = 0.4), size = 7, show.legend = FALSE)+
  guides(color = guide_legend(keyheight = 2.5))

levo_key <- ggplot(levofloxacin_key, aes(x= parameter, y=Odds, color = bacteria, group = bacteria))+
  geom_hline(yintercept = 1, color = "black")+
  geom_point(size = 4)+
  theme_bw()+
  scale_y_continuous(expand= c(0,0), limits = c(0,2.3))+
  xlab("Infection Site")+
  labs(color = "Species")+
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.position = "none")+
  scale_color_viridis_d(option = "D", labels = label_wrap(10))+
  geom_text(
    data = levofloxacin_key,
    aes(x= parameter, y=max(Odds)+0.2, color = bacteria, label = sig),
    position = position_dodgev(height = 0.4), size = 7, show.legend = FALSE)+
  guides(color = guide_legend(keyheight = 2.5))

levo_who <- ggplot(levofloxacin_who, aes(x= parameter, y=Odds, color = bacteria, group = bacteria))+
  geom_hline(yintercept = 1, color = "black")+
  geom_point(size = 4, position = position_dodge(width = 0.1))+
  theme_bw()+
  scale_y_continuous(expand= c(0,0), limits = c(0,5.6))+
  xlab("WHO region")+
  labs(color = "Species")+
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none")+
  scale_color_viridis_d(option = "D", labels = label_wrap(10))+
  geom_text(
    data = levofloxacin_who,
    aes(x= parameter, y=max(Odds)+0.6, color = bacteria, label = sig),
    position = position_dodgev(height = 0.8), size = 7, show.legend = FALSE)+
  guides(color = guide_legend(keyheight = 2.5))

require("patchwork")

combined_plot <- (levo_age | levo_key / levo_who/levo_gen)+plot_layout(guides = "collect",)

ggsave("plots/levofloxacin_mv_results.png", combined_plot, height = 10, width = 19)
