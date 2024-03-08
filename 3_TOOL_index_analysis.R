##### MICAG: Exploration of MIC distribution index comparison
library(tidyverse);library(ggplot2);library(cowplot); library(patchwork)
theme_set(theme_bw())

## Index outputs generated in 1a_TOOL_screening_plots.R
# That code must be run for the chosen characteristic with and without gender 
# for this analysis to have inputs

## What characteristic to look at. (Note: Must match column name)
characteristics <- c("age_group", "key_source") #Options in the default data are: 
#"key_source" # "age_group" # country # income_grp #who_region

for(characteristic in characteristics){
  # read in the data
  index_comparison_gender <-read_csv(paste0("plots/gender_",characteristic,"index_store.csv"))
  index_comparison <- read_csv(paste0("plots/",characteristic,"index_store.csv"))
  
  # Explore data: how many high level and extract max level
  sum_index_gender <- index_comparison_gender %>% 
    group_by(antibiotic, organism_clean, gender, MIC) %>%
    summarise(df_mic = max(dff)) %>% # Get one value per MIC ##### substantial n differences: add , total = sum(Total) here? 
    group_by(antibiotic, organism_clean, gender) %>%
    summarise(mx = max(df_mic), # Max diff for this bug_drug 
              n_big = sum(unique(df_mic) > 0.1)) # Count how many MIC have > 10% diffs
  sum_index <- index_comparison %>%
    group_by(antibiotic, organism_clean, MIC) %>%
    summarise(df_mic = max(dff)) %>% # Get one value per MIC 
    group_by(antibiotic, organism_clean) %>%
    summarise(mx = max(df_mic), # Max diff for this bug_drug 
              n_big = sum(unique(df_mic) > 0.1)) # Count how many MIC have > 10% diffs
  
  # Explore index separately by with / without gender
  g1 <- ggplot(sum_index_gender %>% filter(n_big > 3), aes(y=antibiotic, x = mx, group = interaction(organism_clean, gender))) + geom_point(aes(col = organism_clean, pch = gender)) + 
    ggtitle(characteristic) + 
    scale_x_continuous("Maximum difference in MIC\nacross groupings")
  ggsave(paste0("plots/", characteristic, "index_all.pdf"), height = 10, width = 7)
  
  g2 <- ggplot(sum_index %>% filter(n_big > 3), aes(y=antibiotic, x = mx, group = organism_clean)) + geom_point(aes(col = organism_clean)) + 
    scale_x_continuous("Maximum difference in MIC\nacross groupings")
  
  g1 / g2 + plot_layout(guides = "collect")& theme(legend.position = 'bottom')
  ggsave(paste0("plots/", characteristic, "sep_index.pdf"), height = 10, width = 7)
  
  # Combine with and without gender
  plot_index <- rbind(sum_index_gender, sum_index %>% mutate(gender = "N"))
  
  ggplot(plot_index %>% filter(n_big > 3), aes(x=antibiotic, y = mx, group = interaction(organism_clean, gender))) + 
    geom_point(aes(col = organism_clean, pch = gender), size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
    ggtitle(characteristic) + 
    scale_y_continuous(limits = c(0,0.4), "Maximum difference in MIC\nacross groupings") + 
    scale_shape_discrete("sex",breaks = c("f","m","N"), labels = c("Female","Male","Both")) + 
    geom_hline(yintercept = c(0.1,0.2,0.3), lty = "dashed") + 
    theme(legend.text = element_text(face = "italic"))
  ggsave(paste0("plots/", characteristic, "index_tog.pdf"), height = 10, width = 7)
  
  # Output for use in figures plot 
  write.csv(plot_index, paste0("plots/", characteristic, "index.csv"))
  
  ### Stats for report: what bacteria are in the final differences? 
  high_all <- plot_index %>% filter(n_big > 3)
  
  round(100*table(high_all%>% ungroup() %>% filter(gender == "N") %>% dplyr::select(organism_clean)) / dim(high_all%>% ungroup() %>% filter(gender == "N"))[1],0)
  # In men
  round(100*table(high_all %>% ungroup() %>% filter(gender == "m") %>% dplyr::select(organism_clean)) / dim(high_all %>% filter(gender=="m"))[1],0)
  # In women
  round(100*table(high_all %>% ungroup() %>% filter(gender == "f") %>% dplyr::select(organism_clean)) / dim(high_all %>% filter(gender=="f"))[1],0)
  
  ##### OVER TIME 
  # read in the data
  index_comparison_gender_yr <-read_csv(paste0("plots/year_gender_",characteristic,"index_store.csv"))
  index_comparison_yr <- read_csv(paste0("plots/year_",characteristic,"index_store.csv"))
  
  # Explore data: how many high level and extract max level over time 
  sum_index_gender_yr <- index_comparison_gender_yr %>% 
    group_by(antibiotic, organism_clean, gender, MIC, year) %>%
    summarise(df_mic = max(dff)) %>% # Get one value per MIC 
    group_by(antibiotic, organism_clean, gender, year) %>%
    summarise(mx = max(df_mic), # Max diff for this bug_drug 
              n_big = sum(unique(df_mic) > 0.1)) # Count how many MIC have > 10% diffs
  sum_index_yr <- index_comparison_yr %>%
    group_by(antibiotic, organism_clean, MIC, year) %>%
    summarise(df_mic = max(dff)) %>% # Get one value per MIC 
    group_by(antibiotic, organism_clean, year) %>%
    summarise(mx = max(df_mic), # Max diff for this bug_drug 
              n_big = sum(unique(df_mic) > 0.1)) # Count how many MIC have > 10% diffs
  
  # Exploratory plots: hard to see past data issues and heterogeneity
  ggplot(sum_index_gender_yr %>% filter(n_big > 3), aes(x=antibiotic, y = mx, group = interaction(organism_clean, gender))) + 
    geom_point(aes(col = organism_clean, pch = gender)) + 
    facet_wrap(~year) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
    ggtitle(characteristic) + 
    scale_y_continuous("Maximum difference in MIC\nacross groupings")
  
  ggplot(sum_index_gender_yr %>% filter(n_big > 3), aes(x=antibiotic, y = mx, group = interaction(organism_clean, gender))) + 
    geom_point(aes(col = year, pch = gender)) + 
    facet_wrap(~organism_clean) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
    ggtitle(characteristic) + 
    scale_y_continuous("Maximum difference in MIC\nacross groupings")
  
  ggplot(sum_index_gender_yr, 
         aes(x=antibiotic, y = mx, group = interaction(organism_clean, gender, year))) + 
    geom_line(aes(col = factor(year))) + 
    facet_grid(gender~organism_clean) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
    ggtitle(characteristic) + 
    scale_y_continuous("Maximum difference in MIC\nacross groupings")
  
  # heat map
  gg <- sum_index_gender_yr %>% filter(n_big > 3)
  ggplot(gg, aes(x=year, y = antibiotic, z = mx)) + 
    geom_tile(aes(fill = mx)) + 
    facet_grid(gender~organism_clean) + 
    ggtitle(characteristic) + 
    scale_fill_continuous("Maximum\nindex") + 
    theme(strip.text = element_text(face = "italic"))
  ggsave(paste0("plots/", characteristic, "index_time_heat_map_allbac.pdf"), height = 7, width = 15)
  
  ggplot(gg %>% filter(organism_clean %in% c("Staphylococcus aureus","Escherichia coli")), aes(x=year, y = antibiotic, z = mx)) + 
    geom_tile(aes(fill = mx)) + 
    facet_grid(gender~organism_clean) + 
    ggtitle(characteristic) + 
    scale_fill_continuous("Maximum\nindex") + 
    theme(strip.text = element_text(face = "italic"))
  ggsave(paste0("plots/", characteristic, "index_time_heat_map.pdf"), height = 7, width = 15)
  
}

