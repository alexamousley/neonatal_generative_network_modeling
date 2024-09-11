# Individual generative model analysis and visualizations

# Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
# Email: alexa.mousley@mrc-cbu.cam.ac.uk

# This script runs generalized additive models of the model fit/parameters as
# well as exploring network organization

# 1) Model Fits and Parameters 
#   - Violin plots of best-fit model energies and topological fingerprint dissimilarity
#   - GAMs for energy, eta and gamma
# 2) Global Organization
#   - Linear regressions of simulated x observed global network measures
# 3) Spatial embedding
#   - Linear mixed effects models of local network measures with nodes and participants as random effects

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
library('ggExtra')
library('gratia') 
library('lmerTest')

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
## Observed ##
obs_global_statistics <- readMat('observed_global_statistics.mat')
obs_global_statistics <- as.data.frame(obs_global_statistics$observed.global.statistics)
obs_local_statistics <- readMat('observed_local_statistics.mat')
# Load each local statistic
obs_degree <- as.data.frame(obs_local_statistics$observed.local.statistics[,,1])
obs_betweenness <- as.data.frame(obs_local_statistics$observed.local.statistics[,,2])
obs_clustering <- as.data.frame(obs_local_statistics$observed.local.statistics[,,3])
obs_edge_length <- as.data.frame(obs_local_statistics$observed.local.statistics[,,4])
obs_local_efficiency <- as.data.frame(obs_local_statistics$observed.local.statistics[,,5])
obs_matching <- as.data.frame(obs_local_statistics$observed.local.statistics[,,6])
## Simulated ##
setwd('/set/your/path')                                 # <<<<<<<<< SET
sim_global_statistics <- readMat('individual_simulated_global_statistics.mat')
sim_global_statistics <- as.data.frame(sim_global_statistics$simulated.global.statistics)
sim_local_statistics <- readMat('individual_simulated_local_statistics.mat')
# Load each local statistic
sim_degree <- as.data.frame(sim_local_statistics$simulated.local.statistics[,,1])
sim_betweenness <- as.data.frame(sim_local_statistics$simulated.local.statistics[,,2])
sim_clustering <- as.data.frame(sim_local_statistics$simulated.local.statistics[,,3])
sim_edge_length <- as.data.frame(sim_local_statistics$simulated.local.statistics[,,4])
sim_local_efficiency <- as.data.frame(sim_local_statistics$simulated.local.statistics[,,5])
sim_matching <- as.data.frame(sim_local_statistics$simulated.local.statistics[,,6])

# Load model fit data
energy <- readMat('individual_energy_sorted.mat')
energy <- as.data.frame(energy$Eall.sorted)     # This is all the energies for every run
best_energy <- energy[,1]                       # Take the energy of the best-fit model
params <- readMat('individual_mean_top_params.mat')
params <- as.data.frame(params$mean.top.params)
TFdissimilarity <- readMat('individual_tfdissimilarity.mat')
TFdissimilarity <- as.data.frame(TFdissimilarity$tfdissimilarity)

# Create model dataframe
data <- bind_cols(postconceptional_age,birth_age,sex,head_circum,translation,rotation,best_energy,params)
colnames(data) <- c("postconceptional_age","birth_age","sex","head_circum","translation","rotation",
                    'energy','eta','gamma')

# Create global organization dataframe
data_global <- bind_cols(obs_global_statistics,sim_global_statistics)
colnames(data_global) <- c('obs_modularity','obs_path_length','obs_global_efficiency','obs_rich_connections',
                           'obs_feeder_connections','obs_local_connections','obs_rich_lengths','obs_feeder_lengths',
                           'obs_local_lengths','sim_modularity','sim_path_length','sim_global_efficiency','sim_rich_connections',
                           'sim_feeder_connections','sim_local_connections','sim_rich_lengths','sim_feeder_lengths',
                           'sim_local_lengths')

############################################## (1) Model Fits and Parameters ##############################################

### Violin Plots ###
## Energy
x <- rep(1,nrow(data))                            # Create numbers of participants
energy_df <- bind_cols(x,best_energy)             # Create energy dataframe (participant x best fit energy)
colnames(energy_df) <- c("Participants","Energy") # Name columns
# Plot (Figure 4C)
energy_violin <- energy_df %>%
  ggplot(aes(x=Participants, y = Energy))+
  geom_violin(fill = '#b6e3c4') + 
  theme_classic() +
  stat_summary(fun=mean,  geom="point", size=2, color="black") +
  theme(text=element_text(family="Arial",size=20),axis.title.x = element_blank(),
        axis.text.x = element_blank(),legend.position='none',axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15,colour='black')) 
## TF dissimilarity
TFdiss_df <- bind_cols(x,best_energy)             
colnames(TFdiss_df) <- c("Participants","TFdiss") 
# Plot (Figure 4D)
TF_violin <- TFdiss_df %>%
  ggplot(aes(x=Participants, y = TFdiss))+
  geom_violin(fill = '#b6e3c4') + 
  theme_classic() +
  stat_summary(fun=mean,  geom="point", size=2, color="black") +
  theme(text=element_text(family="Arial",size=20),axis.title.x = element_blank(),
        axis.text.x = element_blank(),legend.position='none',axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15,colour='black'))+labs(y='TF Dissimilarity')

### Energy GAM ###
energy_model <- gam(energy ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                        +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(energy_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(energy_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(energy ~ s(postconceptional_age,bs="cs",k=3)+s(birth_age,bs="cs",k=3)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(energy_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(energy_model)

# Plot energy model (Figure 5A)
reg_plots <- visreg(energy_model,gg=TRUE,type="conditional") 
energy_plot <- ggplot() +
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
  xlim(23, 45.5) + labs(y='Energy',x='Weeks since conception')

### Eta GAM ###
eta_model <- gam(eta ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                               +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(eta_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(eta_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(eta ~ s(postconceptional_age,bs="cs",k=3)+s(birth_age,bs="cs",k=3)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(eta_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(eta_model)

# Plot eta model (Figure 2B)
reg_plots <- visreg(eta_model,gg=TRUE,type="conditional") 
eta_plot <- ggplot() +
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
  xlim(23, 45.5) + labs(y='Eta',x='Weeks since conception')

### Gamma GAM ###
gamma_model <- gam(gamma ~ s(postconceptional_age,bs="cr",k=50)+s(birth_age,bs="cr",k=50)
                              +sex+head_circum+translation+rotation, data=data,method='REML')

# Assess model
gam.check(gamma_model)  # Check 'k' value (if significant, increase k and rerun)
appraise(gamma_model)   # Create plots
# Determine if wiggle is needed (via cs basis function)
model_cs <- gam(gamma ~ s(postconceptional_age,bs="cs",k=3)+s(birth_age,bs="cs",k=3)
                +sex+head_circum+translation+rotation, data=data,method='REML')
p <- par(mfrow=c(2,2))    # Plot and visually inspet if model_cs reduced wiggle
p1 <- plot.gam(gamma_model,rug=T)
p2 <- plot.gam(model_cs,rug=T)

# Print model summary
summary(gamma_model)

# Plot gamma model (Figure 2B)
reg_plots <- visreg(gamma_model,gg=TRUE,type="conditional") 
gamma_plot <- ggplot() +
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
  xlim(23, 45.5) + labs(y='Gamma',x='Weeks since conception')

### Parameter by parameter plots ###

# Create regression plots
ereg_plots <- visreg(eta_model,gg=TRUE,type="conditional")   # Eta
greg_plots <- visreg(gamma_model,gg=TRUE,type="conditional") # Gamma
# Create term indexes
term <- birth_age >= 37
# Group average
terme <- mean(ereg_plots[[2]]$data$y[term]) # Eta
termg <- mean(greg_plots[[2]]$data$y[term]) # Gamma
preterme <- mean(ereg_plots[[2]]$data$y[!term])
pretermg <- mean(greg_plots[[2]]$data$y[!term])

# Plot birth age predicted parameters (Figure 5B)
birth_predicted_plot <- ggplot() +
  geom_point(aes(x = ereg_plots[[2]]$data$y, y = greg_plots[[2]]$data$y, shape = term, colour = term, fill = term), size = 6, alpha = 0.7) +
  geom_point(aes(x = ereg_plots[[2]]$data$y, y = greg_plots[[2]]$data$y, shape = term, colour = term, fill = term), size = 6, alpha = 0.7) +
  geom_point(aes(x = terme, y = termg), colour = 'yellow', shape = 'triangle', size = 8) +
  geom_point(aes(x = preterme, y = pretermg), colour = 'yellow', shape = 'circle', size = 8) +
  scale_fill_manual(values = c("#2d3063", '#02c0d7')) +
  scale_color_manual(values = c("#2d3063", '#02c0d7')) + 
  scale_shape_manual(values = c(21, 24)) +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 50), legend.position = "none",
        axis.title = element_blank()) +
  labs(x = "Eta", y = "Gamma", shape = "Term Group", fill = "Term Group") 
birth_predicted_plot <- ggMarginal(birth_predicted_plot, groupColour = TRUE, groupFill = TRUE, 
                             margins = 'both', type = "density")
# Plot postconceptional age predicted parameters (Figure 5B)
scan_predicted_plot <- ggplot() +
  geom_point(aes(x = ereg_plots[[1]]$data$y, y = greg_plots[[1]]$data$y, colour = ereg_plots[[1]]$data$x, fill = ereg_plots[[1]]$data$x), size = 6, alpha = 0.7, data = reg_plots[[1]]$data) +
  scale_color_gradient(low = '#A44F00', high = "#FFDB1C") +
  theme_classic() +
  theme(text = element_text(family = "Arial", size = 50),
        axis.title = element_blank(),
        legend.position = 'none') +
  labs(x = "Eta", y = "Gamma") 
scan_predicted_plot <- ggMarginal(scan_predicted_plot, groupColour = FALSE, 
                                  margins = 'both', type = "density", fill = "#E69F00")

#################################### (2) Global Organization: Figure 4E #########################################
# Modularity 
modularity_plot <- ggplot(data_global,aes(x=obs_modularity,y=sim_modularity)) + 
  geom_point(alpha=0.5,size=5) +
  scale_color_gradient() +
  geom_smooth(method = 'lm',formula = y ~ x,color='red',size=5) +
  theme_classic() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        axis.line = element_line(linetype='solid',linewidth = 4))
# Global Efficiency
global_efficiency_plot <- ggplot(data_global,aes(x=obs_global_efficiency,y=sim_global_efficiency)) + 
  geom_point(alpha=0.5,size=5) +
  scale_color_gradient() +
  geom_smooth(method = 'lm',formula = y ~ x,color='red',size=5) +
  theme_classic() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        axis.line = element_line(linetype='solid',linewidth = 4))
# Path Length
path_length_plot <- ggplot(data_global,aes(x=obs_path_length,y=sim_path_length)) + 
  geom_point(alpha=0.5,size=5) +
  scale_color_gradient() +
  geom_smooth(method = 'lm',formula = y ~ x,color='red',size=5) +
  theme_classic() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        axis.line = element_line(linetype='solid',linewidth = 4))
# Rich Connections
rich_connections_plot <- ggplot(data_global,aes(x=obs_rich_connections,y=sim_rich_connections)) + 
  geom_point(alpha=0.5,size=5) +
  scale_color_gradient() +
  geom_smooth(method = 'lm',formula = y ~ x,color='red',size=5) +
  theme_classic() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        axis.line = element_line(linetype='solid',linewidth = 4))
# Rich Connection Lengths
rich_lengths_plot <- ggplot(data_global,aes(x=obs_rich_lengths,y=sim_rich_lengths)) + 
  geom_point(alpha=0.5,size=5) +
  scale_color_gradient() +
  geom_smooth(method = 'lm',formula = y ~ x,color='red',size=5) +
  theme_classic() +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        axis.line = element_line(linetype='solid',linewidth = 4))

############################################### (3) Spatial Embedding #################################################

### Degree Linear Model ### 
# Add participant column
obs_degree <- bind_cols(1:nrow(data),obs_degree)
sim_degree <- bind_cols(1:nrow(data),sim_degree)
# Rename columns with node number
colnames(obs_degree) <- c('Participants',1:90)
colnames(sim_degree) <- c('Participants',1:90)
# Pivot dataframes into 3 columns (participants x node x degree)
obs_degree_long <- pivot_longer(obs_degree,col=2:91,
                       names_to=c("node"),values_to = "observed") 
sim_degree_long <- pivot_longer(sim_degree,col=2:91,
                       names_to=c("node"),values_to = "simulated")
# Combine observed and simulated data
degree_all <- bind_cols(obs_degree_long,sim_degree_long[,3])
# Remove subcortical and all 0 gyrus
degree_cortical <- filter(degree_all, node != 71:80)  
# Create participant factor
degree_cortical$ParticipantPlot <- as.factor(degree_cortical$Participants)

# Run
degree_model <- lmer(simulated ~ observed + (1 + observed|Participants) + (1 + observed|node), 
                          data=degree_cortical) # random and slopes intercepts
# Print model summary
summary(degree_model)

# Check model
plot(degree_model, subject.n = NULL)
qqnorm(residuals(degree_model))

# Plot model individual regressions and average (Figure 4F and Supplementary Figure 6)
degree_plot <- degree_cortical %>%
  ggplot(aes(x=observed, y=simulated,shape=ParticipantPlot)) +
  stat_smooth(method='lm',geom='line', alpha=0.1, se=FALSE) +
  geom_smooth(method = 'lm', se = F, aes(group = 1),color='red',size=2)+
  theme_classic()+ theme(axis.text=element_text(size=30),
                         axis.title=element_blank(),
                         legend.position="none",
                         plot.title = element_blank(),
                         text=element_text(family="Arial")) 


### Betweenness Linear Model ###
# Add participant column
obs_betweenness <- bind_cols(1:nrow(data),obs_betweenness)
sim_betweenness <- bind_cols(1:nrow(data),sim_betweenness)
# Rename columns with node number
colnames(obs_betweenness) <- c('Participants',1:90)
colnames(sim_betweenness) <- c('Participants',1:90)
# Pivot dataframes into 3 columns (participants x node x betweenness)
obs_betweenness_long <- pivot_longer(obs_betweenness,col=2:91,
                                    names_to=c("node"),values_to = "observed") 
sim_betweenness_long <- pivot_longer(sim_betweenness,col=2:91,
                                    names_to=c("node"),values_to = "simulated")
# Combine observed and simulated data
betweenness_all <- bind_cols(obs_betweenness_long,sim_betweenness_long[,3])
# Remove subcortical and all 0 gyrus
betweenness_cortical <- filter(betweenness_all, node != 71:80)  
# Create participant factor
betweenness_cortical$ParticipantPlot <- as.factor(betweenness_cortical$Participants)

# Run
betweenness_model <- lmer(simulated ~ observed + (1 + observed|Participants) + (1 + observed|node), 
                         data=betweenness_cortical) # random and slopes intercepts
# Print model summary
summary(betweenness_model)

# Check model
plot(betweenness_model, subject.n = NULL)
qqnorm(residuals(betweenness_model))

# Plot model individual regressions and average (Figure 4F and Supplementary Figure 6)
betweenness_plot <- betweenness_cortical %>%
  ggplot(aes(x=observed, y=simulated,shape=ParticipantPlot)) +
  stat_smooth(method='lm',geom='line', alpha=0.1, se=FALSE) +
  geom_smooth(method = 'lm', se = F, aes(group = 1),color='grey',size=2)+
  theme_classic()+ theme(axis.text=element_text(size=30),
                         axis.title=element_blank(),
                         legend.position="none",
                         plot.title = element_blank(),
                         text=element_text(family="Arial"))

### Clustering Linear Model ###
# Add participant column
obs_clustering <- bind_cols(1:nrow(data),obs_clustering)
sim_clustering <- bind_cols(1:nrow(data),sim_clustering)
# Rename columns with node number
colnames(obs_clustering) <- c('Participants',1:90)
colnames(sim_clustering) <- c('Participants',1:90)
# Pivot dataframes into 3 columns (participants x node x clustering)
obs_clustering_long <- pivot_longer(obs_clustering,col=2:91,
                                    names_to=c("node"),values_to = "observed") 
sim_clustering_long <- pivot_longer(sim_clustering,col=2:91,
                                    names_to=c("node"),values_to = "simulated")
# Combine observed and simulated data
clustering_all <- bind_cols(obs_clustering_long,sim_clustering_long[,3])
# Remove subcortical and all 0 gyrus
clustering_cortical <- filter(clustering_all, node != 71:80)  
# Create participant factor
clustering_cortical$ParticipantPlot <- as.factor(clustering_cortical$Participants)

# Run
clustering_model <- lmer(simulated ~ observed + (1 + observed|Participants) + (1 + observed|node), 
                         data=clustering_cortical) # random and slopes intercepts
# Print model summary
summary(clustering_model)

# Check model
plot(clustering_model, subject.n = NULL)
qqnorm(residuals(clustering_model))

# Plot model individual regressions and average (Figure 4F and Supplementary Figure 6)
clustering_plot <- clustering_cortical %>%
  ggplot(aes(x=observed, y=simulated,shape=ParticipantPlot)) +
  stat_smooth(method='lm',geom='line', alpha=0.1, se=FALSE) +
  geom_smooth(method = 'lm', se = F, aes(group = 1),color='red',size=2)+
  theme_classic()+ theme(axis.text=element_text(size=30),
                         axis.title=element_blank(),
                         legend.position="none",
                         plot.title = element_blank(),
                         text=element_text(family="Arial"))

### Edge Length Linear Model ###
# Add participant column
obs_edge_length <- bind_cols(1:nrow(data),obs_edge_length)
sim_edge_length <- bind_cols(1:nrow(data),sim_edge_length)
# Rename columns with node number
colnames(obs_edge_length) <- c('Participants',1:90)
colnames(sim_edge_length) <- c('Participants',1:90)
# Pivot dataframes into 3 columns (participants x node x edge_length)
obs_edge_length_long <- pivot_longer(obs_edge_length,col=2:91,
                                    names_to=c("node"),values_to = "observed") 
sim_edge_length_long <- pivot_longer(sim_edge_length,col=2:91,
                                    names_to=c("node"),values_to = "simulated")
# Combine observed and simulated data
edge_length_all <- bind_cols(obs_edge_length_long,sim_edge_length_long[,3])
# Remove subcortical and all 0 gyrus
edge_length_cortical <- filter(edge_length_all, node != 71:80)  
# Create participant factor
edge_length_cortical$ParticipantPlot <- as.factor(edge_length_cortical$Participants)

# Run
edge_length_model <- lmer(simulated ~ observed + (1 + observed|Participants) + (1 + observed|node), 
                         data=edge_length_cortical) # random and slopes intercepts
# Print model summary
summary(edge_length_model)

# Check model
plot(edge_length_model, subject.n = NULL)
qqnorm(residuals(edge_length_model))

# Plot model individual regressions and average (Figure 4F and Supplementary Figure 6)
edge_length_plot <- edge_length_cortical %>%
  ggplot(aes(x=observed, y=simulated,shape=ParticipantPlot)) +
  stat_smooth(method='lm',geom='line', alpha=0.1, se=FALSE) +
  geom_smooth(method = 'lm', se = F, aes(group = 1),color='red',size=2)+
  theme_classic()+ theme(axis.text=element_text(size=30),
                         axis.title=element_blank(),
                         legend.position="none",
                         plot.title = element_blank(),
                         text=element_text(family="Arial"))

### Local Efficiency Linear Model ###
# Add participant column
obs_local_efficiency <- bind_cols(1:nrow(data),obs_local_efficiency)
sim_local_efficiency <- bind_cols(1:nrow(data),sim_local_efficiency)
# Rename columns with node number
colnames(obs_local_efficiency) <- c('Participants',1:90)
colnames(sim_local_efficiency) <- c('Participants',1:90)
# Pivot dataframes into 3 columns (participants x node x local_efficiency)
obs_local_efficiency_long <- pivot_longer(obs_local_efficiency,col=2:91,
                                    names_to=c("node"),values_to = "observed") 
sim_local_efficiency_long <- pivot_longer(sim_local_efficiency,col=2:91,
                                    names_to=c("node"),values_to = "simulated")
# Combine observed and simulated data
local_efficiency_all <- bind_cols(obs_local_efficiency_long,sim_local_efficiency_long[,3])
# Remove subcortical and all 0 gyrus
local_efficiency_cortical <- filter(local_efficiency_all, node != 71:80)  
# Create participant factor
local_efficiency_cortical$ParticipantPlot <- as.factor(local_efficiency_cortical$Participants)

# Run
local_efficiency_model <- lmer(simulated ~ observed + (1 + observed|Participants) + (1 + observed|node), 
                         data=local_efficiency_cortical) # random and slopes intercepts
# Print model summary
summary(local_efficiency_model)

# Check model
plot(local_efficiency_model, subject.n = NULL)
qqnorm(residuals(local_efficiency_model))

# Plot model individual regressions and average (Figure 4F and Supplementary Figure 6)
local_efficiency_plot <- local_efficiency_cortical %>%
  ggplot(aes(x=observed, y=simulated,shape=ParticipantPlot)) +
  stat_smooth(method='lm',geom='line', alpha=0.1, se=FALSE) +
  geom_smooth(method = 'lm', se = F, aes(group = 1),color='red',size=2)+
  theme_classic()+ theme(axis.text=element_text(size=30),
                         axis.title=element_blank(),
                         legend.position="none",
                         plot.title = element_blank(),
                         text=element_text(family="Arial"))

### Matching Linear Model ###
# Add participant column
obs_matching <- bind_cols(1:nrow(data),obs_matching)
sim_matching <- bind_cols(1:nrow(data),sim_matching)
# Rename columns with node number
colnames(obs_matching) <- c('Participants',1:90)
colnames(sim_matching) <- c('Participants',1:90)
# Pivot dataframes into 3 columns (participants x node x matching)
obs_matching_long <- pivot_longer(obs_matching,col=2:91,
                                    names_to=c("node"),values_to = "observed") 
sim_matching_long <- pivot_longer(sim_matching,col=2:91,
                                    names_to=c("node"),values_to = "simulated")
# Combine observed and simulated data
matching_all <- bind_cols(obs_matching_long,sim_matching_long[,3])
# Remove subcortical and all 0 gyrus
matching_cortical <- filter(matching_all, node != 71:80)  
# Create participant factor
matching_cortical$ParticipantPlot <- as.factor(matching_cortical$Participants)

# Run
matching_model <- lmer(simulated ~ observed + (1 + observed|Participants) + (1 + observed|node), 
                         data=matching_cortical) # random and slopes intercepts
# Print model summary
summary(matching_model)

# Check model
plot(matching_model, subject.n = NULL)
qqnorm(residuals(matching_model))

# Plot model individual regressions and average (Figure 4F and Supplementary Figure 6)
matching_plot <- matching_cortical %>%
  ggplot(aes(x=observed, y=simulated,shape=ParticipantPlot)) +
  stat_smooth(method='lm',geom='line', alpha=0.1, se=FALSE) +
  geom_smooth(method = 'lm', se = F, aes(group = 1),color='red',size=2)+
  theme_classic()+ theme(axis.text=element_text(size=30),
                         axis.title=element_blank(),
                         legend.position="none",
                         plot.title = element_blank(),
                         text=element_text(family="Arial"))

