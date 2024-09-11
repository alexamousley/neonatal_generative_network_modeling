# Run statistics for propensity-matched analysis

# Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
# Email: alexa.mousley@mrc-cbu.cam.ac.uk

# This script takes the output of the propensity-matched calculate_organization_measures.m code 
# creates box 


# Clear environment
rm(list=ls())

# Load packages
library('tidyverse')
library('mgcv')
library('dplyr')
library('R.matlab')
library('visreg')
library('ggplot2')
library('gridExtra')
library('gratia') 

# Load demographic data
setwd("/imaging/astle/am11/dHCP/published_data/propensity_matched_analysis/demographics/") # <<<<<<<<< SET
birth_age <- readMat('PM_GA_at_birth.mat')
birth_age <- as.data.frame(birth_age$GA.at.birth)
postconceptional_age <- readMat('PM_postconceptional_age.mat')
postconceptional_age <- as.data.frame(postconceptional_age$postconceptional.age)
head_circum <- readMat('PM_head_circumference.mat')
head_circum <- as.data.frame(head_circum$head.circumference)
translation <- readMat('PM_translation.mat')
translation <- as.data.frame(translation$translation)
rotation <- readMat('PM_rotation.mat')
rotation <- as.data.frame(rotation$rotation)
sex <- readMat('PM_sex.mat')
sex <- as.data.frame(sex$sex)

# Load network data
setwd('/imaging/astle/am11/dHCP/published_data/propensity_matched_analysis/derived_data/') # <<<<<<<< SET
global_statistics <- readMat('propensity_global_statistics.mat')
global_statistics <- as.data.frame(global_statistics$propensity.global.statistics)

# Create one dataframe
data <- bind_cols(postconceptional_age,birth_age,sex,head_circum,translation,rotation,global_statistics)
colnames(data) <- c("postconceptional_age","birth_age","sex","head_circum","translation","rotation",
                    'modularity','global_efficiency')
data$Group <- ifelse(data$birth_age < 37 ,'Preterm', 'Term')

#### Modularity ####
# Run t-test
test_result <- t.test(modularity ~ Group, data = data)
# Print the test result
print(test_result)
# Create boxplot
modularity_boxplot <- ggplot(data, aes(x = Group, y = modularity, fill = Group)) +
  geom_boxplot() + theme_classic() +  
  scale_fill_manual(values = c("Preterm" = "#2d3063", "Term" = '#02c0d7')) +
  theme(text = element_text(family = "Arial", size = 60),
        axis.title = element_blank()) + labs(y = "Modularity") + guides(fill = guide_legend(
          title = NULL,
          title.theme = element_blank(),
          label.theme = element_text(size = 20)))
modularity_boxplot

#### Global Efficiency ####
# Run t-test
test_result <- t.test(global_efficiency ~ Group, data = data)
# Print the test result
print(test_result)
# Create a boxplot
ge_boxplot <- ggplot(data, aes(x = Group, y = global_efficiency, fill = Group)) +
  geom_boxplot() + theme_classic() +  
  scale_fill_manual(values = c("Preterm" = "#2d3063", "Term" = '#02c0d7')) +
  theme(text = element_text(family = "Arial", size = 60),
        axis.title = element_blank()) + labs(y = "Global Efficiency") + guides(fill = guide_legend(
          title = NULL,
          title.theme = element_blank(),
          label.theme = element_text(size = 20)))
ge_boxplot
