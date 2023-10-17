##### Functions 


### Plot to give initial exploration plots and csv on index values
plot_generation_MICAG <- function(data, bacteria, groupings, gender_options = T){
  ## data: data table with minimum columns of mic, organism_clean, gender and the user entered characteristics (see below)
  ## bacteria: which organisms to consider for this analysis
  ## groupings: which groupings to analyse MIC differences by
  ## gender_options: baseline include gender differences, can set to F 
  
  ######*********************** RUN ************************#################
  # after specified the items above, just run the whole script and it will 
  # generate the desired plots in a sub-folder called plots, as well as generating values for the index
  
  for(characteristic in groupings){
    for (include_gender in gender_options){
      
      print(paste0("Running for ", characteristic))
      
      # make sure there's a folder to store the plots
      dir.create(file.path("plots"), showWarnings = FALSE)
      
      # empty baseline 
      index_store <- c()
      output_plot <- c()
      
      # Look at patterns in the bacteria with or without gender
      if(include_gender == F){
        for(j in bacteria){
          data_sub <- data[organism_clean == j]
          
          # vector for storing relevant drugs and plots
          drugs <- unique(data_sub$antibiotic)
          
          drugs <- sort(drugs)
          
          plot_store <- list()
          
          #for each of the relevant drugs
          for(i in drugs){
            data_sub_drug <- data_sub[antibiotic == i]
            #count observations by subset
            test <- data_sub_drug[, .N, by = .(mic, get(characteristic) )]
            colnames(test) <- c("MIC", characteristic, "N")
            # also need out of total observations for the age_group/gender
            test2 <- data_sub_drug[, .N, by = .(get(characteristic))]
            colnames(test2) <- c(characteristic, "N")
            # note total number of MIC samples
            tot_samps <- sum(test2$N)
            # combine the two together so can work out proportion
            test[test2, on = c(characteristic), Total := i.N]
            #work out proportion
            test[, prop := N/Total]
            # cumulative sum of proportion (first order)
            test <- test[order(MIC, get(characteristic))]
            for_plot <-test[, cumulative_sum := cumsum(prop), by = c(characteristic)]
            
            # store plot
            for_plot <- for_plot[N>100]
            if(nrow(for_plot)>0){
              temp<- ggplot(for_plot[N>=100], aes(x= MIC, y =cumulative_sum, colour = !!sym(characteristic))) + 
                geom_line()+
                labs(title = paste0("MIC - ", i, paste0(". Tot samples = ", tot_samps)), x = "MIC value", 
                     y = paste0("cumulative proportion of samples by ", characteristic), 
                     colour = characteristic) + 
                scale_x_log10() + 
                theme_linedraw() 
            }
            ### Output 
            output_plot <- rbind(output_plot, for_plot %>% mutate(antibiotic = i, organism_clean = j))
            
            ## Explore index
            if(characteristic == "key_source"){
              for_plot <- for_plot %>% filter(!key_source == "") # remove this from index comparison
            }
            suppressWarnings(index_store <- rbind(index_store, for_plot %>% group_by(MIC) %>% 
                                                    mutate(dff = diff(range(cumulative_sum))) %>% mutate(antibiotic = i, organism_clean = j)))
            # warning when no difference
            plot_store[[i]] <- temp
          }
          
          tiff(paste0("plots/",j , "_", characteristic, "_MICs.tiff"), width = 2500, height = 1500)
          print(cowplot::plot_grid(plotlist =  plot_store) )
          dev.off()
          
        }
        write.csv(index_store, paste0("plots/",characteristic, "index_store.csv"))
        write.csv(output_plot, paste0("plots/",characteristic, "output.csv"))
      }
      
      if(include_gender == T){
        
        for(j in bacteria){
          
          
          data_sub <- data[organism_clean == j]
          # vector for storing relevant drugs and plots
          # vector for storing relevant drugs and plots
          
          drugs <- unique(data_sub$antibiotic)
          drugs <- sort(drugs)
          
          plot_store <- list()
          
          #for each of the relevant drugs
          for(i in drugs){
            data_sub_drug <- data_sub[antibiotic == i]
            #count observations by subset
            test <- data_sub_drug[, .N, by = .(gender, mic, get(characteristic) )]
            colnames(test) <- c("gender","MIC", characteristic, "N")
            # also need out of total observations for the age_group/gender
            test2 <- data_sub_drug[, .N, by = .(gender, get(characteristic))]
            colnames(test2) <- c("gender",characteristic, "N")
            # note total number of MIC samples
            tot_samps <- sum(test2$N)
            # combine the two together so can work out proportion
            test[test2, on = c(characteristic, "gender"), Total := i.N]
            #work out proportion
            test[, prop := N/Total]
            # cumulative sum of proportion (first order)
            test <- test[order(MIC, gender,  get(characteristic))]
            for_plot <-test[, cumulative_sum := cumsum(prop), by = c("gender", characteristic)]
            
            # store plot
            
            if(nrow(for_plot)>0){
              temp<- ggplot(for_plot, aes(x= MIC, y =cumulative_sum, colour = !!sym(characteristic), 
                                          linetype = gender)) + 
                geom_line()+
                labs(title = paste0("MIC by age group - ", i, paste0(". Tot samples = ", tot_samps)), x = "MIC value", 
                     y = paste0("cumulative proportion of samples by ", characteristic), 
                     colour = characteristic) + 
                scale_x_log10() + 
                theme_linedraw() 
            }
            ### Output 
            output_plot <- rbind(output_plot, for_plot %>% mutate(antibiotic = i, organism_clean = j))
            
            ## Explore index
            if(characteristic == "key_source"){
              for_plot <- for_plot %>% filter(!key_source == "") # remove this from index comparison
            }
            suppressWarnings( index_store <- rbind(index_store, for_plot %>% group_by(MIC, gender) %>% 
                                                     mutate(dff = diff(range(cumulative_sum))) %>% mutate(antibiotic = i, organism_clean = j)))
            #Warnings if no difference.
            
            plot_store[[i]] <- temp
          }
          
          tiff(paste0("plots/gender_",j , "_", characteristic, "_MICs.tiff"), width = 2500, height = 1500)
          print(cowplot::plot_grid(plotlist =  plot_store) )
          dev.off()  
        }
        write.csv(index_store, paste0("plots/gender_",characteristic, "index_store.csv"))
        write.csv(output_plot, paste0("plots/gender_",characteristic, "output.csv"))
        
      }
    }
  }
}

### Analysis with time 
plot_generation_bytime_MICAG <- function(data, bacteria, groupings, gender_options = T){
  ## data: data table with minimum columns of mic, organism_clean, gender and the user entered characteristics (see below)
  ## bacteria: which organisms to consider for this analysis
  ## groupings: which groupings to analyse MIC differences by
  ## gender_options: baseline include gender differences, can set to F 
  
  ######*********************** RUN ************************#################
  # after specified the two above items, can just run the whole script and it will 
  # generate the desired plots
  for(characteristic in characteristics){
    
    print(paste0("Running for ", characteristic))
    
    for (include_gender in include_gender_options){
      
      # make sure there's a folder to store the plots
      dir.create(file.path("plots"), showWarnings = FALSE)
      index_store <- c()
      output_plot <- c()
      
      # Look at patterns in three bacteria with or without gender
      if(include_gender == F){
        for(j in bacteria_to_use){
          data_sub <- full_data[organism == j]
          
          # vector for storing relevant drugs and plots
          drugs <- unique(data_sub$antibiotic)
          
          drugs <- sort(drugs)
          
          plot_store <- list()
          
          #for each of the relevant drugs
          for(i in drugs){
            data_sub_drug <- data_sub[antibiotic == i]
            #count observations by subset
            test <- data_sub_drug[, .N, by = .(mic, year, get(characteristic) )]
            colnames(test) <- c("MIC","year", characteristic, "N")
            # also need out of total observations for the age_group/gender
            test2 <- data_sub_drug[, .N, by = .(get(characteristic), year)]
            colnames(test2) <- c(characteristic,"year", "N")
            # note total number of MIC samples
            tot_samps <- sum(test2$N)
            # combine the two together so can work out proportion
            test[test2, on = c(characteristic, "year"), Total := i.N]
            #work out proportion
            test[, prop := N/Total]
            # cumulative sum of proportion (first order)
            test <- test[order(MIC, year, get(characteristic))]
            for_plot <-test[, cumulative_sum := cumsum(prop), by = c("year",characteristic)]
            
            # store plot
            for_plot <- for_plot[N>100]
            if(nrow(for_plot)>0){
              temp<- ggplot(for_plot[N>=100], aes(x= MIC, y =cumulative_sum, colour = !!sym(characteristic), group = year)) + 
                geom_line()+
                labs(title = paste0("MIC - ", i, paste0(". Tot samples = ", tot_samps)), x = "MIC value", 
                     y = paste0("cumulative proportion of samples by ", characteristic), 
                     colour = characteristic) + 
                scale_x_log10() + 
                theme_linedraw() 
            }
            ### Output 
            output_plot <- rbind(output_plot, for_plot %>% mutate(antibiotic = i, organism = j))
            
            ## Explore index
            if(characteristic == "key_source"){
              for_plot <- for_plot %>% filter(!key_source == "") # remove this from index comparison
            }
            suppressWarnings( index_store <- rbind(index_store, for_plot %>% group_by(MIC) %>%
                                                     mutate(dff = diff(range(cumulative_sum))) %>% mutate(antibiotic = i, organism = j)))
            #warnings if no difference
            
            plot_store[[i]] <- temp
          }
          
          tiff(paste0("plots/year_",j , "_", characteristic, "_MICs.tiff"), width = 2500, height = 1500)
          print(cowplot::plot_grid(plotlist =  plot_store) )
          dev.off()  
          
          
        }
        write.csv(index_store, paste0("plots/year_",characteristic, "index_store.csv"))
        write.csv(output_plot, paste0("plots/year_",characteristic, "output.csv"))
      }
      
      if(include_gender == T){
        for(j in bacteria_to_use){
          
          data_sub <- full_data[organism == j]
          # vector for storing relevant drugs and plots
          # vector for storing relevant drugs and plots
          drugs <- unique(data_sub$antibiotic)
          drugs <- sort(drugs)
          
          plot_store <- list()
          
          #for each of the relevant drugs
          for(i in drugs){
            data_sub_drug <- data_sub[antibiotic == i]
            #count observations by subset
            test <- data_sub_drug[, .N, by = .(gender, mic, year, get(characteristic) )]
            colnames(test) <- c("gender","MIC","year", characteristic, "N")
            # also need out of total observations for the age_group/gender
            test2 <- data_sub_drug[, .N, by = .(gender, year, get(characteristic))]
            colnames(test2) <- c("gender","year",characteristic, "N")
            # note total number of MIC samples
            tot_samps <- sum(test2$N)
            # combine the two together so can work out proportion
            test[test2, on = c(characteristic, "gender","year"), Total := i.N]
            #work out proportion
            test[, prop := N/Total]
            # cumulative sum of proportion (first order)
            test <- test[order(MIC, gender, year, get(characteristic))]
            for_plot <-test[, cumulative_sum := cumsum(prop), by = c("gender","year", characteristic)]
            
            # store plot
            
            if(nrow(for_plot)>0){
              temp<- ggplot(for_plot, aes(x= MIC, y =cumulative_sum, colour = !!sym(characteristic), 
                                          linetype = gender)) + 
                geom_line()+
                labs(title = paste0("MIC by age group - ", i, paste0(". Tot samples = ", tot_samps)), x = "MIC value", 
                     y = paste0("cumulative proportion of samples by ", characteristic), 
                     colour = characteristic) + 
                scale_x_log10() + 
                theme_linedraw() 
            }
            ### Output 
            output_plot <- rbind(output_plot, for_plot %>% mutate(antibiotic = i, organism = j))
            
            ## Explore index
            if(characteristic == "key_source"){
              for_plot <- for_plot %>% filter(!key_source == "") # remove this from index comparison
            }
            
            suppressWarnings(index_store <- rbind(index_store, for_plot %>% group_by(MIC,year, gender) %>% 
                                                    mutate(dff = diff(range(cumulative_sum))) %>% mutate(antibiotic = i, organism = j)))
            # warning when no difference
            
            plot_store[[i]] <- temp
          }
          
          tiff(paste0("plots/year_gender_",j , "_", characteristic, "_MICs.tiff"), width = 2500, height = 1500)
          print(cowplot::plot_grid(plotlist =  plot_store) )
          dev.off()  
          
        }
        write.csv(index_store, paste0("plots/year_gender_",characteristic, "index_store.csv"))
        write.csv(output_plot, paste0("plots/year_gender_",characteristic, "output.csv"))
      }
    }
  }
}