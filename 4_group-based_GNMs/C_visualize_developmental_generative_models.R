# Developmental generative model analysis and visualizations

# Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
# Email: alexa.mousley@mrc-cbu.cam.ac.uk

# This script runs generalized additive models on the connection data from the developmental generative models

# 1) Connection Counts
#   - Rich, feeder, local connections counts plotted separately
# 2) Average connection lengths
# 3) Connection lengths for each connection type
#   - Rich, feeder, local connection lengths plotted separately

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
setwd('/set/your/path/')                       # <<<<<< SET
term_counts <- readMat('mean_term_counts.mat')
term_counts <- as.data.frame(term_counts$mean.term.counts)
preterm_counts <- readMat('mean_preterm_counts.mat')
preterm_counts <- as.data.frame(preterm_counts$mean.preterm.counts)
term_lengths <- readMat('total_term_lengths.mat')
term_lengths <- as.data.frame(t(term_lengths$total.term.lengths))
preterm_lengths <- readMat('total_preterm_lengths.mat')
preterm_lengths <- as.data.frame(t(preterm_lengths$total.preterm.lengths))
term_lengths_per_type <- readMat('mean_term_lengths.mat')
term_lengths_per_type <- as.data.frame(term_lengths_per_type$mean.term.lengths)
preterm_lengths_per_type <- readMat('mean_preterm_lengths.mat')
preterm_lengths_per_type <- as.data.frame(preterm_lengths_per_type$mean.preterm.lengths)

# Create iteration and group variables
iteration <- as.data.frame(rep(1:nrow(preterm_counts)))
group <- c(rep(1,361),rep(2,361))

# Create dataframes
models <- bind_cols(rbind(iteration,iteration),group,rbind(term_counts,preterm_counts),
                    rbind(term_lengths,preterm_lengths),
                    rbind(term_lengths_per_type,preterm_lengths_per_type))
colnames(models) <- c('iteration','group','rich_counts','feeder_counts','local_counts',
                      'lengths','rich_lengths','feeder_lengths','local_lengths')

grouped_models <- bind_cols(iteration,term_counts,preterm_counts,term_lengths,preterm_lengths,
                    term_lengths_per_type,preterm_lengths_per_type)
colnames(grouped_models) <- c('iteration','term_rich_count','term_feeder_count','term_local_count',
                           'preterm_rich_count','preterm_feeder_count','preterm_local_count',
                           'term_length','preterm_length',
                           'term_rich_length','term_feeder_length','term_local_length',
                           'preterm_rich_length','preterm_feeder_length','preterm_local_length')

############################################## (1) Connection Counts #####################################################

### Rich Connections ###
# Run initial model
rich_counts_model <- gam(rich_counts ~ s(iteration,bs="cr",k=10)+group,
                        data=models,method='REML')
# Print model summary
summary(rich_counts_model) # If 'group' is significant, make seperate models for visualization

# Make seperate models
term_rich_counts_model <- gam(term_rich_count ~ s(iteration,bs="cr",k=10),
                              data=grouped_models,method='REML')
term_reg_plots <- visreg(term_rich_counts_model,gg=TRUE,type="conditional") 
preterm_rich_counts_model <- gam(preterm_rich_count ~ s(iteration,bs="cr",k=10),
                              data=grouped_models,method='REML')
preterm_reg_plots <- visreg(preterm_rich_counts_model,gg=TRUE,type="conditional") 

# Plot (Figure 6B)
rich_connections_plot <- ggplot()+
  geom_line(aes(x=term_reg_plots$layers[[3]]$data$x, y=term_reg_plots$layers[[3]]$data$y),colour='red',linewidth=4)+
  geom_line(aes(x=preterm_reg_plots$layers[[3]]$data$x, y=preterm_reg_plots$layers[[3]]$data$y),colour="red",linetype=3,linewidth=4)+
  theme_classic() + 
  theme(text=element_text(family="Arial",size=60),axis.title=element_blank())

### Feeder Connections ###
# Run initial model
feeder_counts_model <- gam(feeder_counts ~ s(iteration,bs="cr",k=10)+group,
                         data=models,method='REML')
# Print model summary
summary(feeder_counts_model) # If 'group' is significant, make seperate models for visualization

# Make seperate models
term_feeder_counts_model <- gam(term_feeder_count ~ s(iteration,bs="cr",k=10),
                              data=grouped_models,method='REML')
term_reg_plots <- visreg(term_feeder_counts_model,gg=TRUE,type="conditional") 
preterm_feeder_counts_model <- gam(preterm_feeder_count ~ s(iteration,bs="cr",k=10),
                                 data=grouped_models,method='REML')
preterm_reg_plots <- visreg(preterm_feeder_counts_model,gg=TRUE,type="conditional") 

# Plot (Figure 6B)
feeder_connections_plot <- ggplot()+
  geom_line(aes(x=term_reg_plots$layers[[3]]$data$x, y=term_reg_plots$layers[[3]]$data$y),colour='#CCCC00',linewidth=4)+
  geom_line(aes(x=preterm_reg_plots$layers[[3]]$data$x, y=preterm_reg_plots$layers[[3]]$data$y),colour="#CCCC00",linetype=3,linewidth=4)+
  theme_classic() + 
  theme(text=element_text(family="Arial",size=60),axis.title=element_blank())

### Local Connections ###
# Run initial model
local_counts_model <- gam(local_counts ~ s(iteration,bs="cr",k=10)+group,
                         data=models,method='REML')
# Print model summary
summary(local_counts_model) # If 'group' is significant, make seperate models for visualization

# Make seperate models
term_local_counts_model <- gam(term_local_count ~ s(iteration,bs="cr",k=10),
                              data=grouped_models,method='REML')
term_reg_plots <- visreg(term_local_counts_model,gg=TRUE,type="conditional") 
preterm_local_counts_model <- gam(preterm_local_count ~ s(iteration,bs="cr",k=10),
                                 data=grouped_models,method='REML')
preterm_reg_plots <- visreg(preterm_local_counts_model,gg=TRUE,type="conditional") 

# Plot (Figure 6B)
local_connections_plot <- ggplot()+
  geom_line(aes(x=term_reg_plots$layers[[3]]$data$x, y=term_reg_plots$layers[[3]]$data$y),colour='black',linewidth=4)+
  geom_line(aes(x=preterm_reg_plots$layers[[3]]$data$x, y=preterm_reg_plots$layers[[3]]$data$y),colour="black",linetype=3,linewidth=4)+
  theme_classic() + 
  theme(text=element_text(family="Arial",size=60),axis.title=element_blank())

############################################## (2) Total Connection Lengths ###############################################

# Run initial model
lengths_model <- gam(lengths ~ s(iteration,bs="cr",k=10)+group,
                         data=models,method='REML')
# Print model summary
summary(lengths_model) # If 'group' is significant, make seperate models for visualization

# Make seperate models
term_lengths_model <- gam(term_length ~ s(iteration,bs="cr",k=10),
                              data=grouped_models,method='REML')
term_reg_plots <- visreg(term_lengths_model,gg=TRUE,type="conditional") 
preterm_lengths_model <- gam(preterm_length ~ s(iteration,bs="cr",k=10),
                                 data=grouped_models,method='REML')
preterm_reg_plots <- visreg(preterm_lengths_model,gg=TRUE,type="conditional") 

# Plot (Figure 6C)
length_plot <- ggplot()+
  geom_line(aes(x=term_reg_plots$layers[[3]]$data$x, y=term_reg_plots$layers[[3]]$data$y),colour='#02c0d7',linewidth=4)+
  geom_line(aes(x=preterm_reg_plots$layers[[3]]$data$x, y=preterm_reg_plots$layers[[3]]$data$y),colour="#2d3063",linewidth=4)+
  theme_classic() + 
  theme(text=element_text(family="Arial",size=60),axis.title=element_blank())

############################################## (3) Connection Lengths Per Type ############################################

### Rich Connections ###
# Run initial model
rich_lengths_model <- gam(rich_lengths ~ s(iteration,bs="cr",k=10)+group,
                         data=models,method='REML')
# Print model summary
summary(rich_lengths_model) # If 'group' is significant, make seperate models for visualization

# Make seperate models
term_rich_lengths_model <- gam(term_rich_length ~ s(iteration,bs="cr",k=10),
                              data=grouped_models,method='REML')
term_rich_plots <- visreg(term_rich_lengths_model,gg=TRUE,type="conditional") 
preterm_rich_lengths_model <- gam(preterm_rich_length ~ s(iteration,bs="cr",k=10),
                                 data=grouped_models,method='REML')
preterm_rich_plots <- visreg(preterm_rich_lengths_model,gg=TRUE,type="conditional") 

### Feeder Connections ###
# Run initial model
feeder_lengths_model <- gam(feeder_lengths ~ s(iteration,bs="cr",k=10)+group,
                           data=models,method='REML')
# Print model summary
summary(feeder_lengths_model) # If 'group' is significant, make seperate models for visualization

# Make seperate models
term_feeder_lengths_model <- gam(term_feeder_length ~ s(iteration,bs="cr",k=10),
                                data=grouped_models,method='REML')
term_feeder_plots <- visreg(term_feeder_lengths_model,gg=TRUE,type="conditional") 
preterm_feeder_lengths_model <- gam(preterm_feeder_length ~ s(iteration,bs="cr",k=10),
                                   data=grouped_models,method='REML')
preterm_feeder_plots <- visreg(preterm_feeder_lengths_model,gg=TRUE,type="conditional") 

### Local Connections ###
# Run initial model
local_lengths_model <- gam(local_lengths ~ s(iteration,bs="cr",k=10)+group,
                          data=models,method='REML')
# Print model summary
summary(local_lengths_model) # If 'group' is significant, make seperate models for visualization

# Make seperate models
term_local_lengths_model <- gam(term_local_length ~ s(iteration,bs="cr",k=10),
                               data=grouped_models,method='REML')
term_local_plots <- visreg(term_local_lengths_model,gg=TRUE,type="conditional") 
preterm_local_lengths_model <- gam(preterm_local_length ~ s(iteration,bs="cr",k=10),
                                  data=grouped_models,method='REML')
preterm_local_plots <- visreg(preterm_local_lengths_model,gg=TRUE,type="conditional") 

# Plot (Figure 6D)
length_by_connection_plot <- ggplot()+
  # Rich
  geom_line(aes(x=term_rich_plots$layers[[3]]$data$x, y=term_rich_plots$layers[[3]]$data$y),colour='red',linewidth=4)+
  geom_line(aes(x=preterm_rich_plots$layers[[3]]$data$x, y=preterm_rich_plots$layers[[3]]$data$y),colour="red",linetype=3,linewidth=4)+
  # Feeder
  geom_line(aes(x=term_feeder_plots$layers[[3]]$data$x, y=term_feeder_plots$layers[[3]]$data$y),colour='#CCCC00',linewidth=4)+
  geom_line(aes(x=preterm_feeder_plots$layers[[3]]$data$x, y=preterm_feeder_plots$layers[[3]]$data$y),colour='#CCCC00',linetype=3,linewidth=4)+
  # Local
  geom_line(aes(x=term_local_plots$layers[[3]]$data$x, y=term_local_plots$layers[[3]]$data$y),colour='black',linewidth=4)+
  geom_line(aes(x=preterm_local_plots$layers[[3]]$data$x, y=preterm_local_plots$layers[[3]]$data$y),colour="black",linetype=3,linewidth=4)+
  theme_classic() + 
  theme(text=element_text(family="Arial",size=60),axis.title=element_blank())



