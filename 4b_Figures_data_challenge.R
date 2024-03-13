# MICAG: This script creates the figures used in the report as part of the Data Challenge

library(data.table);library(ggplot2);library(cowplot); library(patchwork); library(tidyverse); library(viridis)
theme_set(theme_bw(base_size = 16))
# read in the data
#full_data <- as.data.table(read.csv("data/full_data.csv"))
#source("overlapping_drugs.R")


# Which characteristics explored? 
characteristics <- c("age_group", "key_source") # "age_group"

############################# Get data  ######
output_data <- c()
output_index <- c()

for(i in characteristics){
  print(i)
  
  ###### Output no time 
  output1 <- read.csv(paste0("plots/gender_",i, "output.csv")) %>% mutate(charac = i) %>% rename("charac_value" = i)
  output2 <- read.csv(paste0("plots/",i, "output.csv")) %>% mutate(gender = "N", charac = i) %>% rename("charac_value" = i)
  output_data <- rbind(rbind(output_data, 
                             output1 %>% dplyr::select("gender", "MIC", "charac", "charac_value", "N", "Total", "prop", "cumulative_sum", 
                                                "antibiotic", "organism_clean")), 
                       output2 %>% dplyr::select("gender", "MIC", "charac", "charac_value", "N", "Total", "prop", "cumulative_sum", 
                                          "antibiotic", "organism_clean"))
  
  ###### Index comparison
  index_comparison <- read.csv(paste0("plots/", i, "index.csv")) %>% mutate(charac = i)
  output_index <- rbind(output_index, index_comparison)
  
}

### EXAMPLES
combinations <- as.data.frame(rbind(c("Staphylococcus aureus", "levofloxacin"),
                                    c("Escherichia coli", "levofloxacin"),
                                    c("Staphylococcus aureus", "meropenem"),
                                    c("Escherichia coli", "meropenem"),
                                    c("Staphylococcus aureus", "ampicillin"),
                                    c("Escherichia coli", "ampicillin")))
colnames(combinations) <- c("organism_clean","antibiotic")


# Grab just the data wanted for the examples
plot_data <- c()
plot_index <- c()

for(i in 1:nrow(combinations)){
  plot_data <- rbind(plot_data,
                     output_data %>% filter(organism_clean == combinations[i,"organism_clean"],
                                            antibiotic == combinations[i,"antibiotic"]))
  
  plot_index <- rbind(plot_index, 
                      output_index %>% filter(organism_clean == combinations[i,"organism_clean"],
                                              antibiotic == combinations[i,"antibiotic"])) 
}

## make names shorter
wsa <- which(plot_data$organism_clean == "Staphylococcus aureus")
plot_data[wsa,"organism_clean"] <- "S. aureus"
wec <- which(plot_data$organism_clean == "Escherichia coli")
plot_data[wec,"organism_clean"] <- "E. coli"
unique(plot_data$organism_clean)

#################################### Age only  ######
plot_age <- plot_data %>% filter(charac == "age_group")
plot_age$charac_value <- factor(plot_age$charac_value, 
                                levels = c("0 to 2 Years","3 to 12 Years", "13 to 18 Years",
                                           "19 to 64 Years", "65 to 84 Years", "85 and Over"))


ggplot(plot_age %>% filter(gender == "N"), 
       aes(x=MIC, y = cumulative_sum, group = charac_value)) + 
  geom_line(aes(col = charac_value)) + 
  facet_grid(organism_clean ~ antibiotic, scales = "free") + 
  scale_x_log10("MIC", labels = scales::comma) + 
  scale_y_continuous("Cumulative proportion of\nisolates tested") + 
  scale_color_discrete("Age group")
ggsave("plots/age_only.pdf")

# Age and gender 
g1 <- ggplot(plot_age %>% filter(!gender == "N"), 
             aes(x=MIC, y = cumulative_sum, group = interaction(gender,charac_value))) + 
  geom_line(aes(col = charac_value, lty = gender)) + 
  facet_grid(organism_clean ~ antibiotic, scales = "free") + 
  scale_x_log10("MIC", labels = scales::comma) + 
  ggtitle("Age and Sex") + 
  scale_y_continuous("Cumulative proportion of\nisolates tested") + 
  #scale_color_discrete("Age group") + 
  scale_linetype_discrete("Sex", labels = c("Female","Male"), breaks = c("f","m")) + 
  theme(legend.position = "bottom", strip.text = element_text(face = "italic")) + 
  scale_colour_viridis("Age group", option="D", discrete = TRUE)
ggsave("plots/age_sex.pdf")

# output_index <- data.table(output_index)
# output_index[antibiotic == "trimethoprim sulfa", antibiotic := "trimethroprim\nsulfa"] 
# output_index[antibiotic == "piperacillin tazobactam", antibiotic := "piperacillin\ntazobactam"]
# output_index[antibiotic == "meropenem vaborbactam", antibiotic := "meropenem\nvaborbactam"]
# output_index[antibiotic == "cefoperazone sulbactam", antibiotic := "piperacillin\ntazobactam"]
# output_index[antibiotic == "ceftolozane tazobactam", antibiotic := "cefoperazone\nsulbactam"]


g2 <- ggplot(output_index %>% filter(charac == "age_group",n_big > 3),
             aes(y=antibiotic, x = mx, group = interaction(organism_clean, gender))) + 
  geom_point(aes(col = organism_clean, pch = gender), size = 3) + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  ggtitle("Age and Sex") + 
  #scale_color_discrete("Organism") + 
  scale_x_continuous(limits = c(0.1,0.4),"Maximum difference in MIC across groupings") + 
  scale_shape_discrete("Sex",breaks = c("f","m","N"), labels = c("Female","Male","Both")) + 
  geom_vline(xintercept = c(0.1,0.2,0.3), lty = "dashed") + 
  theme(legend.text = element_text(face = "italic")) + 
  scale_colour_viridis("Organism", option="D", discrete = TRUE)
ggsave("plots/index_age_sex.pdf")

######################## Source  ######
ggplot(plot_data %>% filter(charac == "key_source", gender == "N"), 
       aes(x=MIC, y = cumulative_sum, group = charac_value)) + 
  geom_line(aes(col = charac_value)) + 
  facet_grid(organism_clean ~ antibiotic, scales = "free") + 
  scale_x_log10("MIC") + 
  scale_y_continuous("Cumulative proportion of\nisolates tested") + 
  scale_color_discrete("Age group")
ggsave("plots/source.pdf")

# Source and gender 
g3 <- ggplot(plot_data %>% filter(charac == "key_source", !charac_value == "", !gender == "N"), 
             aes(x=MIC, y = cumulative_sum, group = interaction(gender,charac_value))) + 
  geom_line(aes(col = charac_value, lty = gender)) + 
  facet_grid(organism_clean ~ antibiotic, scales = "free") + 
  scale_x_log10("MIC", labels = scales::comma) +
  ggtitle("Infection site") + 
  guides(lty = "none") +
  scale_y_continuous("Cumulative proportion of\nisolates tested") + 
 # scale_color_discrete("Infection site") + 
  scale_linetype_discrete("Sex", labels = c("Female","Male"), breaks = c("f","m")) + 
  theme(legend.position = "bottom", strip.text = element_text(face = "italic")) + 
  scale_colour_viridis("Infection site", option="D", discrete = TRUE)
ggsave("plots/source_sex.pdf")

g4 <- ggplot(output_index %>% filter(charac == "key_source",n_big > 3), 
             aes(y=antibiotic, x = mx, group = interaction(organism_clean, gender))) + 
  geom_point(aes(col = organism_clean, pch = gender), size = 3) + 
 # theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  ggtitle("Infection site") + 
 # scale_color_discrete("Organism") + 
  scale_x_continuous(limits = c(0.1,0.3), "Maximum difference in MIC across groupings") + 
  scale_shape_discrete("Sex",breaks = c("f","m","N"), labels = c("Female","Male","Both")) + 
  geom_vline(xintercept = c(0.1,0.2,0.3), lty = "dashed") + 
  theme(legend.text = element_text(face = "italic")) + 
  scale_colour_viridis("Organism", option="D", discrete = TRUE)
ggsave("plots/index_key_source.pdf")

#### Index figure
a <- g2 + g4 +  theme(legend.position = "none")
b <- g1 + g3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
b / a + plot_layout(heights = c(1, 1.8)) + plot_annotation(tag_levels = 'A')
#ggsave("plots/fig2.pdf", width = 20, height = 11)

legend1 <- get_legend(g2 + theme(legend.position = "bottom"))
legend2 <- get_legend(g1 + theme(legend.position = "bottom", 
                                 legend.justification = "right"))
legend3 <- get_legend(g3 + theme(legend.position = "bottom", legend.justification = "left", 
                                 legend.direction = "horizontal"))


top_plot <- plot_grid(g1 + theme(legend.position = "none") , 
                      g3 + theme(legend.position = "none"), 
                      labels = c("A", "B"), 
                      ncol = 2)
middle_plot <- plot_grid(legend, labels = "")
bottom_plot <- plot_grid(g2 + theme(legend.position = "none") , 
                         g4 + theme(legend.position = "none"), 
                          labels = c("C", "D"), 
                         ncol = 2 )
very_bottom <- plot_grid(legend2, legend3, rel_widths = c(0.6,0.5))
together <- plot_grid(top_plot, middle_plot, bottom_plot, very_bottom, ncol =1,
                      rel_heights = c(1,0.1,1.4, 0.2)) + 
  theme(plot.background = element_rect(fill = "white", colour = NA))
print(together)
ggsave("plots/fig2.tiff", width = 20, height = 11)

### Figure for index graphic
ggplot(plot_age %>% filter(gender == "N", antibiotic == "levofloxacin", organism_clean == "Staphylococcus aureus"), 
             aes(x=MIC, y = cumulative_sum, group = interaction(gender,charac_value))) + 
  geom_point(aes(col = charac_value), size = 10, pch = "x") + 
  geom_line(aes(col = charac_value), size = 2) + 
  facet_grid(organism_clean ~ antibiotic, scales = "free") + 
  scale_x_log10("MIC", labels = scales::comma) + 
  scale_y_continuous("Cumulative proportion of\nisolates tested") + 
  scale_color_discrete("Age group") + 
  scale_linetype_discrete("Sex", labels = c("Female","Male"), breaks = c("f","m")) + 
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text.x = element_blank(), strip.text.y = element_blank())




##### OVER TIME ######
# read in the data
### Get data
output_datat <- c()

for(i in characteristics){
  print(i)
  
  ###### Output with time 
  output1t <- read.csv(paste0("plots/year_gender_",i, "output.csv")) %>% mutate(charac = i) %>% rename("charac_value" = i)
  output2t <- read.csv(paste0("plots/year_",i, "output.csv")) %>% mutate(gender = "N", charac = i) %>% rename("charac_value" = i)
  output_datat <- rbind(rbind(output_datat, 
                              output1t %>% dplyr::select("year","gender", "MIC", "charac", "charac_value", "N", "Total", "prop", "cumulative_sum", 
                                                  "antibiotic", "organism_clean")), 
                        output2t %>% dplyr::select("year","gender", "MIC", "charac", "charac_value", "N", "Total", "prop", "cumulative_sum", 
                                            "antibiotic", "organism_clean"))
}

### Bug/drug combos
combinations <- as.data.frame(rbind(c("Staphylococcus aureus", "levofloxacin"),
                                    c("Escherichia coli", "levofloxacin"),
                                    c("Staphylococcus aureus", "meropenem"),
                                    c("Escherichia coli", "meropenem"),
                                    c("Staphylococcus aureus", "ampicillin"),
                                    c("Escherichia coli", "ampicillin")))
colnames(combinations) <- c("organism_clean","antibiotic")


# Grab just the data wanted for the examples
plot_datat <- c()

for(i in 1:nrow(combinations)){
  plot_datat <- rbind(plot_datat,
                      output_datat %>% filter(organism_clean == combinations[i,"organism_clean"],
                                              antibiotic == combinations[i,"antibiotic"]))
}

# Index 
### Over time: for MICAG analysis most interested in age_group comparisons with time
characteristic <- "age_group"
index_comparison_gender_yr <-read_csv(paste0("plots/year_gender_",characteristic,"index_store.csv"))

# Explore data: how many high level and extract max level over time 
sum_index_gender_yr <- index_comparison_gender_yr %>% 
  group_by(antibiotic, organism_clean, gender, MIC, year) %>%
  summarise(df_mic = max(dff)) %>% # Get one value per MIC 
  group_by(antibiotic, organism_clean, gender, year) %>%
  summarise(mx = max(df_mic), # Max diff for this bug_drug 
            n_big = sum(unique(df_mic) > 0.1)) # Count how many MIC have > 10% diffs

###### heat maps of index over time 
# Plot with pre-selected bug-drug combos
gg <- sum_index_gender_yr %>% filter(n_big > 3)
ggplot(gg, aes(x=year, y = antibiotic, z = mx)) + 
  geom_tile(aes(fill = mx)) + 
  facet_grid(gender~organism_clean) + 
  ggtitle(characteristic) + 
  scale_fill_continuous("Maximum\nindex") + 
  theme(strip.text = element_text(face = "italic"))
ggsave(paste0("plots/", characteristic, "index_time_heat_map_allbac.pdf"), height = 7, width = 15)

# Plot with only only 2 big-drug combos
gg$antibiotic <- tools::toTitleCase(as.character(gg$antibiotic)) 
gg$antibiotic <- factor(gg$antibiotic, levels=rev(unique(gg$antibiotic)))

g1t <- ggplot(gg %>% filter(organism_clean %in% c("Staphylococcus aureus","Escherichia coli")), aes(x=year, y = antibiotic, z = mx)) + 
  geom_tile(aes(fill = mx)) + 
  facet_grid(gender~organism_clean) + 
  ggtitle("A") + 
  labs(x = "Year", y = "Antibiotics")+
  scale_fill_viridis("Maximum\nindex", option = "D", ) +
  theme(strip.text = element_text(face = "italic", size = 10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10))
ggsave(paste0("plots/", characteristic, "index_time_heat_map.pdf"), height = 10, width = 10)



###### isolates per time plot
samples_store <- data.table(read.csv(paste0("plots/year_gender_",characteristic, "samples_store.csv")))

samples_store$age_group <- factor(samples_store$age_group, 
                                  levels = c("0 to 2 Years","3 to 12 Years", "13 to 18 Years",
                                             "19 to 64 Years", "65 to 84 Years", "85 and Over"))
title_to_use <- expression(paste("B - ", italic("S. aureus"), "- levofloxacin"))

g_samples <- ggplot(samples_store, aes(x=Year, y=N, fill = age_group))+
  geom_col()+
  labs(x = "Year",
       y = "Number of samples",
       fill = "Age group", 
       title = title_to_use) +
  scale_fill_viridis(discrete = TRUE, option = "H") +
  facet_grid(MIC~Gender)  + 
  theme(strip.text = element_text(face = "italic", size = 8),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10))

g1t + g_samples
ggsave(paste0("plots/", characteristic, "time_figure3.pdf"), height = 10, width = 20)


###### Line plots over time ######
plot_datat_staphlevo <- plot_datat %>% filter(antibiotic == "levofloxacin",
                                              organism_clean == "Staphylococcus aureus", 
                                              charac == characteristic) %>%
  ungroup %>% 
  group_by(gender, MIC, charac_value) %>% mutate(total_iso = sum(Total))

g2t <- ggplot(plot_datat_staphlevo %>% filter(total_iso > 1000), 
       aes(x = year, y = cumulative_sum, colour = charac_value, 
           group = interaction(gender, charac_value))) + 
  geom_line(aes(linetype = factor(gender))) + 
  facet_wrap(~MIC) + 
  scale_x_log10()
ggsave(paste0("plots/", characteristic, "time_cumulative.pdf"), height = 7, width = 15)

ggplot(plot_datat_staphlevo %>% filter(total_iso > 10000), 
       aes(x = year, y = cumulative_sum, colour = charac_value, 
           group = interaction(gender, charac_value))) + 
  geom_line(aes(linetype = factor(gender))) + 
  facet_wrap(~MIC) + 
  scale_x_log10()

g1t + g2t & theme(legend.position = "bottom")
ggsave(paste0("plots/", characteristic, "time_figure.pdf"), height = 7, width = 15)



###### OPTIONAL: time plot for N and prop for specific bug-drug combos ######

full_data <- as.data.table(read.csv("data/full_data.csv"))
bacteria_to_use <- unique(full_data$organism_clean) 

### combos
combinations <- as.data.frame(rbind(c("Staphylococcus aureus", "levofloxacin"),
                                    c("Escherichia coli", "levofloxacin"),
                                    c("Staphylococcus aureus", "meropenem"),
                                    c("Escherichia coli", "meropenem"),
                                    c("Staphylococcus aureus", "ampicillin"),
                                    c("Escherichia coli", "ampicillin")))
colnames(combinations) <- c("organism_clean","antibiotic")

# Grab just the data wanted for the examples
plot_datat <- c()

for(i in 1:nrow(combinations)){
  plot_datat <- rbind(plot_datat,
                      full_data %>% filter(organism_clean == combinations[i,"organism_clean"],
                                           antibiotic == combinations[i,"antibiotic"]))
                              }

#for loop to produce one graph per bug-drug combo
for(j in bacteria_to_use){
  data_sub <- plot_datat[organism_clean == j]
  
  # vector for storing relevant drugs and plots
  drugs <- unique(data_sub$antibiotic)
  
  drugs <- sort(drugs)
  
  plot_store <- list()
  
  for(i in drugs){
    data_sub_drug <- data_sub[antibiotic == i]
    test <- data_sub_drug[, .N, by = .(mic,year,age_group )]
    colnames(test) <- c("MIC", "year","age_group", "N")
    # also need out of total observations for the age_group/gender
    test2 <- data_sub_drug[, .N, by = .((mic))]
    colnames(test2) <- c("MIC", "N")
    # note total number of MIC samples
    tot_samps <- sum(test2$N)
    # combine the two together so can work out proportion
    test <- test[test2, on = c("MIC"), Total := i.N]
    #work out proportion
    for_plottt <- test[, prop := N/Total]
    
    if(nrow(for_plottt)>0){
    b1 <- ggplot(for_plottt, aes(x=as.factor(year),y=N,fill=as.factor(MIC))) +
      geom_col(position = "identity") 
    
    b2 <- ggplot(for_plottt, aes(x=as.factor(year),y=prop,color=as.factor(MIC),group=as.factor(MIC)))+
      geom_line()  
   
    b3 <- ggplot(for_plottt, aes(x=as.factor(year),y=N,fill=age_group)) +
      geom_col(position = "identity") +
      facet_wrap(~MIC)
    
    b4 <- ggplot(for_plottt, aes(x=as.factor(year),y= prop,color=age_group)) +
      geom_line(aes(group=interaction(age_group,MIC))) +
      facet_wrap(~MIC) 
    
    fullplot <- ((b1/b2) | (b3/b4))
    
    tiff(paste0("plots/",j , "_", i, "_MICs.tiff"), width = 2500, height = 1500)
    print(fullplot, nrow = 4,ncol=3, guides = 'collect')
    dev.off()
    }

  }
  
}
 

