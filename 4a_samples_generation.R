# Vivli plots for MIC by characteristics/groupings over time (and gender) (stacked bars ver.)
library(data.table);library(ggplot2);library(cowplot); library(dplyr)

# Read in the data
full_data <- as.data.table(read.csv("data/full_data.csv"))
# Ensure there's a folder to store the plots
dir.create(file.path("plots"), showWarnings = FALSE)
# ---------------------------------------------------------------------

### CHOOSE BUG-DRUG COMBINATION
bug <- "Staphylococcus aureus"
drug <- "levofloxacin"
### SPECIFY WHICH CHARACTERISTIC/GROUPING TO STUDY (must match column name)
characteristic <- "age_group"
### INCLUDE GENDER? (T for yes, F for no)
include_gender <- T
### CHOOSE THRESHOLD NUMBER OF YEARLY SAMPLES FOR INCLUDING MEASURES
thresh <- 0

# ---------------------------------------------------------------------

# Set up store
plot_store <- list()
samples_store <- c()

# Filter data by bug-drug and lift all reported MIC values
data_sub <- full_data[organism==bug]
data_sub_drug <- data_sub[antibiotic==drug]
miclist <- sort(unique(data_sub_drug$mic))

if(include_gender==F){
  # Loop to plot time graph for each mic
  for(x in 1:length(miclist)){
    # Filter by MIC
    data_sub_drug_m <- data_sub_drug[mic==miclist[x]]
    # Count samples with that MIC by characteristic/grouping AND year
    test_m_1 <- data_sub_drug_m[,.N,by=.(get(characteristic),year)]
    colnames(test_m_1) <- c(characteristic,"Year","N")
    # Include data for plotting if more samples in year than 'threshold'
    for_plot_m <- test_m_1[N>thresh]
    # Plot proportions of samples with that MIC per year by characteristic/grouping
    if(nrow(for_plot_m)>0){
      temp_m <- ggplot(for_plot_m, aes(x=Year, y=N, fill = !!sym(characteristic)))+
        geom_col(color="black")+
        labs(title = paste0("MIC = ", miclist[x],"mg/L"), x = "Year",
             y = "Number of samples",
             fill = characteristic)+
        scale_x_continuous(breaks=2004:2022,limits=c(2003,2023))+
        scale_y_continuous()+
        theme_linedraw()
      
      samples_store <- rbind(samples_store, for_plot_m  %>% mutate(antibiotic = drug, organism = bug, MIC = paste0(miclist[x],"mg/L")))
      
      # Store the plot for that MIC
      plot_store[[x]] <- temp_m
    }
  }
  
  # Combine plots for each MIC into final figure output
  tiff(paste0("plots/",bug,"_",drug,"_",characteristic,"_MICs_over_time_stacked_ver.tiff"),width=2500,height=1500)
  print(cowplot::plot_grid(plotlist=plot_store))
  dev.off()
  write.csv(samples_store, paste0("plots/year_gender_",characteristic, "samples_store.csv"))
  
}

if(include_gender==T){
  # Loop to plot time graph for each mic
  for(x in 1:length(miclist)){
    # Filter by MIC
    data_sub_drug_m <- data_sub_drug[mic==miclist[x]]
    # Count samples with that MIC by characteristic/grouping AND year
    test_m_1 <- data_sub_drug_m[,.N,by=.(gender,get(characteristic),year)]
    colnames(test_m_1) <- c("Gender",characteristic,"Year","N")
    # Include data for plotting if more samples in year and gender than 'threshold'
    for_plot_m <- test_m_1[N>thresh]
    # Plot proportions of samples with that MIC per year by characteristic/grouping
    if(nrow(for_plot_m)>0){
      temp_m <- ggplot(for_plot_m, aes(x=Year, y=N, fill = !!sym(characteristic)))+
        geom_col(color="black")+
        labs(title = paste0("MIC = ", miclist[x],"mg/L"), x = "Year",
             y = "Number of samples",
             fill = characteristic)+
        scale_x_continuous(breaks=2004:2022,limits=c(2003,2023))+
        scale_y_continuous()+
        facet_wrap(vars(Gender),nrow=1)+
        theme_linedraw()+
        theme(axis.text.x = element_text(angle=90)) + 
        lims(y = c(0,3200))
      
      samples_store <- rbind(samples_store, for_plot_m %>% mutate(antibiotic = drug, organism = bug, MIC = paste0(miclist[x],"mg/L")))
      
      # Store the plot for that MIC
      plot_store[[x]] <- temp_m
    }
  }
  
  # Combine plots for each MIC into final figure output
  tiff(paste0("plots/",bug,"_",drug,"_",characteristic,"_MICs_over_time_stacked_ver_by_gender.tiff"),width=2500,height=1500)
  print(cowplot::plot_grid(plotlist=plot_store))
  dev.off()
  write.csv(samples_store, paste0("plots/year_gender_",characteristic, "samples_store.csv"))
  
}

