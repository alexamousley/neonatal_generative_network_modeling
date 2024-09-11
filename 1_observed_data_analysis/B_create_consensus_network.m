%% Create consensus network
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script takes binarized networks and calculates a representative 
consensus network. With parameters as-is, the resulting 'consensus' variable 
should be exactly equal to the published 'consensus_network'.

The consensus network is calcuated by a function found here: 
https://www.brainnetworkslab.com/coderesources

The function is documented in this publication:
Betzel, R. F., Griffa, A., Hagmann, P., & Mišić, B. (2019). 
Distance-dependent consensus thresholds for generating group-representative
structural brain networks. Network neuroscience, 3(2), 475-496.

%}

%% Add paths and load data
clear; clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Load binarized networks
load('binarized_connectomes.mat'); % or load('density_controlled_binarized_connectomes.mat') for density-controlled analysis
% Load region labels
load('regions.mat');
% Load Euclidean distance
load('euclidean_distance.mat');

%% Create consensus network

% Reshape network (region x region x participants)
A = shiftdim(binarized_connectomes, 4);

% Identify hemispheres
hemiid = zeros(size(regions));           % Initialize variable for hemisphere ID
L_index = find(contains(regions, "L"));  % Find nodes in the left hemisphere
R_index = find(contains(regions, "R"));  % Find nodes in the right hemisphere
hemiid(L_index) = 1;                     % Set left hemisphere nodes to '1'
hemiid(R_index) = 2;                     % Set right hemisphere nodes to '2'

nbins = 40;                              % Choose the number of distance bins

% Create consensus network
[G, Gc] = fcn_group_bins(A, edistance, hemiid, nbins);

consensus = G;  % The distance-based consensus network was selected for further analysis

