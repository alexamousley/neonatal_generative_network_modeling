%% Calculate organizational measure for propensity-matched sample
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script takes binarized networks and calculates global
organizational measures using the Brain Connectivity Toolbox
(https://sites.google.com/site/bctnet/home?authuser=0)

Toolbox publication:
Rubinov, M., & Sporns, O. (2010). Complex network measures of brain 
connectivity: uses and interpretations. Neuroimage, 52(3), 1059-1069.

%}

%% Add paths and load data
clear;clc;

% Add paths
run('/imaging/astle/am11/dHCP/neonatal_generative_modelling/set_paths.m');

% Load binarized networks
load('propensity_binarized_connectomes.mat'); 
% Load atlas euclidean distance
load('euclidean_distance.mat');

nsub = size(binarized_connectomes,1); % Define number of participants


%% Calculate organizational measures

% Define names of chosen local and global measures as strings
globalmeasures = string({'Modularity','Global efficiency'});

% Initialize arrays
observed_global_statistics = zeros(nsub, 2);   
            
% Loop over participants
for sub = 1:nsub
    % Extract the observed network for the current participant
    observed_network = squeeze(binarized_connectomes(sub,:,:));
    
    %%% Global statistics %%%
    % Maximum modularity
    [~, observed_global_statistics(sub, 1)] = modularity_und(observed_network);
    % Characteristic path length and global efficiency
    dist = distance_bin(observed_network);
    [~, observed_global_statistics(sub, 2)] = charpath(dist, 0, 0);
end
