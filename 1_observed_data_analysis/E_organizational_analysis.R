# Run generalized additive models

# Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
# Email: alexa.mousley@mrc-cbu.cam.ac.uk

# This script takes the output of the calculate_organization_measures.m code for OBSERVED network measures 
# and runs generalized additive models. This code works for both the observed average density networks and 
# the density-controlled networks. This code also includes the code to run the term-equivalent age 
# analysis (removing all infants scanned before term-equivalent age).

# 1) Time until scan
# 2) Density GAM
# 3) Global Orgnaization GAMs
#   - Modularity
#   - Global Efficiency
#   - Characteristic Path Length
#   - Rich Connections
#   - Rich Connection Lengths
#   - Rich Connections Proportional to Total Connections

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
setwd('/set/your/path')                                 # <<<<<<<<< SET
birth_age <- readMat('GA_at_birth.mat')
birth_age <- as.data.frame(birth_age$GA.at.birth)
postconceptional_age <- readMat('postconceptional_age.mat')
postconceptional_age <- as.data.frame(postconceptional_age$postconceptional.age)
head_circum <- readMat('head_circumference.mat')
head_circum <- as.data.frame(head_circum$head.circumference)
translation <- readMat('translation.mat')
translation <- as.data.frame(translation$translation)
rotation <- readMat('rotation.mat')
rotation <- as.data.frame(rotation$rotation)
sex <- readMat('sex.mat')
sex <- as.data.frame(sex$sex)

# Load network data
setwd('/set/your/path')                                 # <<<<<<<<< SET
density <- readMat('observed_density.mat')    # Average density
#density <- readMat('observed_density10.mat') # Density-controlled
density <- as.data.frame(density$observed.density)
global_statistics <- readMat('observed_global_statistics.mat')    # Average density
#global_statistics <- readMat('observed_global_statistics10.mat') # Density-controlled
global_statistics <- as.data.frame(global_statistics$observed.global.statistics)

# Load each local statistic
local_statistics <- readMat('observed_local_statistics.mat')    # Average density
#local_statistics <- readMat('observed_local_statistics10.mat') # Density-controlled
degree <- as.data.frame(local_statistics$observed.local.statistics[,,1])
betweenness <- as.data.frame(local_statistics$observed.local.statistics[,,2])
clustering <- as.data.frame(local_statistics$observed.local.statistics[,,3])
edge_length <- as.data.frame(local_statistics$observed.local.statistics[,,4])
local_efficiency <- as.data.frame(local_statistics$observed.local.statistics[,,5])
matching <- as.data.frame(local_statistics$observed.local.statistics[,,6])

# Create one dataframe
data <- bind_cols(postconceptional_age,birth_age,sex,head_circum,translation,rotation,density,global_statistics)
colnames(data) <- c("postconceptional_age","birth_age","sex","head_circum","translation","rotation","density",
                        'modularity','characteristic_path_length','global_efficiency','rich_connections',
                    'feeder_connections','local_connections','rich_lengths','feeder_lengths','local_lengths')

## For term-equivalent age analysis, remove all infants scanned earlier that 37 GA
# data <- data[data$postconceptional_age >= 37, ]
# data$Group <- ifelse(data_matched$birth_age < 37, 'Preterm', 'Term')
## NOTE: If you choose to run the term-equivalent age analysis, remove the PCA regression on plots to create
## the graphs in Supplementary Figure 3.

############################################## (1) Time until scan ##############################################
# Calculate the difference between time of birth and time of scan 
difference <- data$postconceptional_age-data$birth_age

# Plot (Figure 1B)
birth_scan_corr <- ggplot(data,aes(x=birth_age,y=difference)) + 
  geom_point(alpha=0.5,size=5) +
  scale_color_gradient() +
  geom_smooth(method = 'lm',formula = y ~ x,color='red',size=2) +
  theme_classic() + xlim(min(birth_age),max(postconceptional_age))+
  theme(text=element_text(family="Arial",size=60),axis.title=element_blank()) +
  geom_vline(xintercept = 37, linetype="longdash", 
             color = "grey", size=2)

############################################## (2) Density GAM #################################################

# Run GAM
density_model <- gam(density ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                     +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(density_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(density_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(density ~ s(postconceptional_age,bs="cs",k=50)+s(birth_age,bs="cs",k=50)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(density_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(density_model)

# Plot density model (Figure 1C)
reg_plots <- visreg(density_model,gg=TRUE,type="conditional") 
density_plot <- ggplot() +
  # Plot points, regression line and confidence intervals of postconceptional age
  geom_polygon(fill = "gray85", data = reg_plots[[1]]$layers[[1]]$data, aes(x = x, y = y * 100), alpha = 0.7) +
  geom_line(data = reg_plots[[1]]$layers[[3]]$data, aes(x = x, y = y * 100), colour = "#E69F00", linewidth = 2) +
  geom_point(data = reg_plots[[1]]$data, aes(x = x, y = y * 100), colour = "#E69F00", alpha = 0.5) +
  # Plot points, regression line and confidence intervals of birth age
  geom_polygon(fill = "gray85", data = reg_plots[[2]]$layers[[1]]$data, aes(x = x, y = y * 100), alpha = 0.7) +
  geom_line(data = reg_plots[[2]]$layers[[3]]$data, aes(x = x, y = y * 100), colour = "#0072B2", linewidth = 2) +
  geom_point(data = reg_plots[[2]]$data, aes(x = x, y = y * 100), colour = "#0072B2", alpha = 0.5) +
  scale_color_manual(name = 'Age Measure',
                     breaks = c('Postconceptional Age', 'Birth age'),
                     values = c('Postconceptional Age' = '#E69F00', 'Birth Age' = '#0072B2')) +
  theme_classic() + 
  theme(text = element_text(family = "Arial", size = 60)) +
  geom_vline(xintercept = 37, linetype = "longdash", color = "grey", size = 2) +
  xlim(23, 45.5) + labs(y='Density (%)',x='Weeks since conception')

############################################## (3) Global Organization GAMs ##############################################
### Modularity GAM ###
modularity_model <- gam(modularity ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                     +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(modularity_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(modularity_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(modularity ~ s(postconceptional_age,bs="cs",k=50)+s(birth_age,bs="cs",k=50)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(modularity_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(modularity_model)

# Plot modularity model (Figure 2A)
reg_plots <- visreg(modularity_model,gg=TRUE,type="conditional") 
modularity_plot <- ggplot() +
  # Plot points, regression line and confidence intervals of postconceptional age
  geom_polygon(fill = "gray85", data = reg_plots[[1]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[1]]$layers[[3]]$data, aes(x = x, y = y), colour = "#E69F00", linewidth = 2) +
  geom_point(data = reg_plots[[1]]$data, aes(x = x, y = y), colour = "#E69F00", alpha = 0.5) +
  # Plot points, regression line and confidence intervals of birth age
  geom_polygon(fill = "gray85", data = reg_plots[[2]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[2]]$layers[[3]]$data, aes(x = x, y = y), colour = "#0072B2", linewidth = 2) +
  geom_point(data = reg_plots[[2]]$data, aes(x = x, y = y), colour = "#0072B2", alpha = 0.5) +
  scale_color_manual(name = 'Age Measure',
                     breaks = c('Postconceptional Age', 'Birth age'),
                     values = c('Postconceptional Age' = '#E69F00', 'Birth Age' = '#0072B2')) +
  theme_classic() + 
  theme(text = element_text(family = "Arial", size = 60)) +
  geom_vline(xintercept = 37, linetype = "longdash", color = "grey", size = 2) +
  xlim(23, 45.5) + labs(y='Modularity',x='Weeks since conception')

### Global efficiency GAM ###
global_efficiency_model <- gam(global_efficiency ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                        +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(global_efficiency_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(global_efficiency_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(global_efficiency ~ s(postconceptional_age,bs="cs",k=3)+s(birth_age,bs="cs",k=3)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(global_efficiency_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(global_efficiency_model)

# Plot global_efficiency model (Figure 2B)
reg_plots <- visreg(global_efficiency_model,gg=TRUE,type="conditional") 
global_efficiency_plot <- ggplot() +
  # Plot points, regression line and confidence intervals of postconceptional age
  geom_polygon(fill = "gray85", data = reg_plots[[1]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[1]]$layers[[3]]$data, aes(x = x, y = y), colour = "#E69F00", linewidth = 2) +
  geom_point(data = reg_plots[[1]]$data, aes(x = x, y = y), colour = "#E69F00", alpha = 0.5) +
  # Plot points, regression line and confidence intervals of birth age
  geom_polygon(fill = "gray85", data = reg_plots[[2]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[2]]$layers[[3]]$data, aes(x = x, y = y), colour = "#0072B2", linewidth = 2) +
  geom_point(data = reg_plots[[2]]$data, aes(x = x, y = y), colour = "#0072B2", alpha = 0.5) +
  scale_color_manual(name = 'Age Measure',
                     breaks = c('Postconceptional Age', 'Birth age'),
                     values = c('Postconceptional Age' = '#E69F00', 'Birth Age' = '#0072B2')) +
  theme_classic() + 
  theme(text = element_text(family = "Arial", size = 60)) +
  geom_vline(xintercept = 37, linetype = "longdash", color = "grey", size = 2) +
  xlim(23, 45.5) + labs(y='Global Efficiency',x='Weeks since conception')

### Path Length GAM ###
path_length_model <- gam(characteristic_path_length ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                               +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(path_length_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(path_length_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(characteristic_path_length ~ s(postconceptional_age,bs="cs",k=3)+s(birth_age,bs="cs",k=3)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(path_length_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(path_length_model)

# Plot path_length model (Figure 2B)
reg_plots <- visreg(path_length_model,gg=TRUE,type="conditional") 
path_length_plot <- ggplot() +
  # Plot points, regression line and confidence intervals of postconceptional age
  geom_polygon(fill = "gray85", data = reg_plots[[1]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[1]]$layers[[3]]$data, aes(x = x, y = y), colour = "#E69F00", linewidth = 2) +
  geom_point(data = reg_plots[[1]]$data, aes(x = x, y = y), colour = "#E69F00", alpha = 0.5) +
  # Plot points, regression line and confidence intervals of birth age
  geom_polygon(fill = "gray85", data = reg_plots[[2]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[2]]$layers[[3]]$data, aes(x = x, y = y), colour = "#0072B2", linewidth = 2) +
  geom_point(data = reg_plots[[2]]$data, aes(x = x, y = y), colour = "#0072B2", alpha = 0.5) +
  scale_color_manual(name = 'Age Measure',
                     breaks = c('Postconceptional Age', 'Birth age'),
                     values = c('Postconceptional Age' = '#E69F00', 'Birth Age' = '#0072B2')) +
  theme_classic() + 
  theme(text = element_text(family = "Arial", size = 60)) +
  geom_vline(xintercept = 37, linetype = "longdash", color = "grey", size = 2) +
  xlim(23, 45.5) + labs(y='Global Efficiency',x='Weeks since conception')

### Rich connection counts GAM ###
rich_connections_model <- gam(rich_connections ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                               +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(rich_connections_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(rich_connections_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(rich_connections ~ s(postconceptional_age,bs="cs",k=3)+s(birth_age,bs="cs",k=3)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(rich_connections_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(rich_connections_model)

# Plot rich_connections model (Figure 2B)
reg_plots <- visreg(rich_connections_model,gg=TRUE,type="conditional") 
rich_connections_plot <- ggplot() +
  # Plot points, regression line and confidence intervals of postconceptional age
  geom_polygon(fill = "gray85", data = reg_plots[[1]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[1]]$layers[[3]]$data, aes(x = x, y = y), colour = "#E69F00", linewidth = 2) +
  geom_point(data = reg_plots[[1]]$data, aes(x = x, y = y), colour = "#E69F00", alpha = 0.5) +
  # Plot points, regression line and confidence intervals of birth age
  geom_polygon(fill = "gray85", data = reg_plots[[2]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[2]]$layers[[3]]$data, aes(x = x, y = y), colour = "#0072B2", linewidth = 2) +
  geom_point(data = reg_plots[[2]]$data, aes(x = x, y = y), colour = "#0072B2", alpha = 0.5) +
  scale_color_manual(name = 'Age Measure',
                     breaks = c('Postconceptional Age', 'Birth age'),
                     values = c('Postconceptional Age' = '#E69F00', 'Birth Age' = '#0072B2')) +
  theme_classic() + 
  theme(text = element_text(family = "Arial", size = 60)) +
  geom_vline(xintercept = 37, linetype = "longdash", color = "grey", size = 2) +
  xlim(23, 45.5) + labs(y='Connections (#)',x='Weeks since conception')

### Rich connection lengths GAM ###
rich_lengths_model <- gam(rich_lengths ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                              +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(rich_lengths_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(rich_lengths_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(rich_lengths ~ s(postconceptional_age,bs="cs",k=3)+s(birth_age,bs="cs",k=3)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(rich_lengths_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(rich_lengths_model)

# Plot rich_lengths model (Figure 2B)
reg_plots <- visreg(rich_lengths_model,gg=TRUE,type="conditional") 
rich_lengths_plot <- ggplot() +
  # Plot points, regression line and confidence intervals of postconceptional age
  geom_polygon(fill = "gray85", data = reg_plots[[1]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[1]]$layers[[3]]$data, aes(x = x, y = y), colour = "#E69F00", linewidth = 2) +
  geom_point(data = reg_plots[[1]]$data, aes(x = x, y = y), colour = "#E69F00", alpha = 0.5) +
  # Plot points, regression line and confidence intervals of birth age
  geom_polygon(fill = "gray85", data = reg_plots[[2]]$layers[[1]]$data, aes(x = x, y = y), alpha = 0.7) +
  geom_line(data = reg_plots[[2]]$layers[[3]]$data, aes(x = x, y = y), colour = "#0072B2", linewidth = 2) +
  geom_point(data = reg_plots[[2]]$data, aes(x = x, y = y), colour = "#0072B2", alpha = 0.5) +
  scale_color_manual(name = 'Age Measure',
                     breaks = c('Postconceptional Age', 'Birth age'),
                     values = c('Postconceptional Age' = '#E69F00', 'Birth Age' = '#0072B2')) +
  theme_classic() + 
  theme(text = element_text(family = "Arial", size = 60)) +
  geom_vline(xintercept = 37, linetype = "longdash", color = "grey", size = 2) +
  xlim(23, 45.5) + labs(y='Average Lengths',x='Weeks since conception')

### Rich connection lengths proportional to total counts GAM ###
# Initialize
total_connections <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
rich_proportion <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
# Iterate over each row
for (sub in 1:nrow(data)) {
  # Calculate the total connections
  total_connections[sub, 1] <- sum(data[sub, 11:13])
  
  # Calculate the rich proportion
  rich_proportion[sub, 1] <- data[sub, 14] / total_connections[sub, 1]
}

# Create new dataframe
data_prop <- bind_cols(postconceptional_age,birth_age,sex,head_circum,translation,rotation,rich_proportion)
colnames(data_prop) <- c("postconceptional_age","birth_age","sex","head_circum","translation","rotation",'rich_proportion')

# Run model
rich_prop_model <- gam(rich_proportion ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                          +sex+head_circum+translation+rotation, data=data_prop,method='REML')

# Assess model
gam.check(rich_prop_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(rich_prop_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(rich_proportion ~ s(postconceptional_age,bs="cs",k=3)+s(birth_age,bs="cs",k=3)
                +sex+head_circum+translation+rotation, data=data_prop,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(rich_prop_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(rich_prop_model)

