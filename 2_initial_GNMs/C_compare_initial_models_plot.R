# Create model comparision plot (Figure 4A)

# Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
# Email: alexa.mousley@mrc-cbu.cam.ac.uk

# This script takes the best 10 models (based on the lowest 10 energies) for the 13 different generative models

# Clear environment
rm(list=ls())

# Load packages
library('R.matlab')
library('ggplot2')
library('reshape')

# Load data
setwd('/set/your/path/')                             # <<<<<<< SET
energies <- readMat('consensus_Eall_sorted.mat')
energies <- as.data.frame(t(energies$consensus.Eall.sorted))   # This variable will be runs x model

# Define model colors and names
modelcolors <- c("#4b86b4","#b6e3c4","#86becb","#aa98a9")
modelabbs <- c("Spatial","H-Neigh","H-Match","C-Avg","C-Min","C-Max","C-Diff","C-Prod",
               "D-Avg","D-Min","D-Max","D-Diff","D-Prod")
modelcats <- c(1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4)

# Reshape energy
energies_melt <- melt(energies)
energies_melt$modeltype <- as.factor(rep(modelcats, each = 10))

energy_comparision_plot <- energies_melt %>%
  ggplot(aes(x = variable, y = value, color = modeltype)) +
  geom_point(size=6) +
  scale_color_manual(labels = c("Spatial", "Homophily","Clustering","Degree"),values = modelcolors)+
  scale_x_discrete(labels = modelabbs)+
  theme_classic() +
  theme(legend.position = c(0.99, 0.99),    
        legend.justification = c(1, 1),legend.box.background = element_rect(color = "black", fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10,face='bold',hjust = 0.5),text = element_text(family = "Arial"),
        axis.text.y = element_text(size=30,colour='black'),axis.text.x=element_text(size=10,colour='black'),
        axis.title = element_blank()) +
  labs(colour = "Model Category") +guides(colour = guide_legend(nrow = 2, ncol = 2, byrow = TRUE))

