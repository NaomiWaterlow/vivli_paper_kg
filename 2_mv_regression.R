##### MICAG regression analysis for specified drugs
# Specify drug / bug combination to run regression analysis
library(data.table); library(MASS)

# read in data (can use alternative)
full_data <- as.data.table(read.csv("data/full_data.csv"))

######*********************** SPECIFY ************************#################
# specify which bacteria and antibiotic of interest
# NOTE - must match spelling in the data table
target_antibiotic <- "levofloxacin"
target_bug <- 'Staphylococcus aureus'

# Specify which age group to make base age group (by putting it first)
# Default is 19-64 years (adults)
full_data$age_group <- factor(full_data$age_group, levels = c(
  "19 to 64 Years", 
  "0 to 2 Years",
  "3 to 12 Years", 
  "13 to 18 Years",
  "65 to 84 Years", 
  "85 and Over"
))

######*********************** RUN ************************#################
# after specified the two above items, can just run the whole script and it will 
# generate the regression and save a table with the coefficient values

# subset to just look at one bug-dryg
sub_data <- full_data[antibiotic == target_antibiotic & 
                         organism_clean == target_bug]

# convert to categorical for each value
sub_data[, mic_cat_all := round(log(mic)/log(2))]
add_to_make_0 <- -min(unique(sub_data$mic_cat_all))
sub_data[, mic_cat_all := factor(add_to_make_0 + round(log(mic)/log(2)))]

sub_data[, year_scaled := year - 2004]
# try running a proportional odds ordinal model
ord_mod <- polr(mic_cat_all ~ age_group + gender  + key_source + year_scaled , data = sub_data, 
                Hess = T)
summary(ord_mod)

summary_table <- coef(summary(ord_mod))
pval <- pnorm(abs(summary_table[, "t value"]),lower.tail = FALSE)* 2
summary_table <- cbind(summary_table, "p value" = round(pval,3))
summary_table <- as.data.frame(summary_table)
summary_table$parameter <- rownames(summary_table)
summary_table <- data.table(summary_table)
summary_table[, Odds := exp(Value)]

summary_table <- summary_table[1:12,c("parameter", "Value", "Std. Error", "p value", "Odds")]

summary_table$Value <- round(summary_table$Value, 3)
summary_table$`Std. Error` <- round(summary_table$`Std. Error`, 3)
summary_table$Odds <- round(summary_table$Odds, 3)
write.csv(summary_table, file = paste0("regresssion_coefficients_", target_antibiotic, 
                                       "_", target_bug,".csv"))
