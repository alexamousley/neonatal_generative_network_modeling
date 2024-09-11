%% Identify rich club nodes and connection types
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script takes the consensus network and identifies rich club nodes. 
Then, it takes binarized networks and identifies the type of connections 
based on the pre-defined rich club nodes. As the parameters are set, the 
resulting 'rich_club_nodes' variable should be exactly equal to the 
published 'rich_club_nodes' variable.

%}

%% Add paths and load data
clear; clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Load consensus network
load('consensus_network.mat'); % or load('density_controlled_consensus_network.mat') for density-controlled analysis
% Load binarized networks
load('binarized_connectomes.mat'); % or load('density_controlled_binarized_connectomes.mat') for density-controlled analysis
% Load atlas Euclidean distance
load('euclidean_distance.mat');

null_sample = 1000;  % Choose the number of null models

%% Calculate observed rich club coefficients
max_degree = max(degrees_und(consensus_network));       % Identify the maximum degree of the consensus network
consensus_rich_coefs = rich_club_bu(consensus_network); % Calculate rich club coefficients of the consensus network

%% Generate null models
% Initialize
null_models = zeros(null_sample, 90, 90);
null_rich_coefs = [];

% Loop through number of null samples
for i = 1:null_sample
    null = randmio_und(consensus_network, 50); % Create null model
    null_models(i, :, :) = null;                       % Save null model
    null_rich_coefs(i, :) = rich_club_bu(null);        % Calculate rich club coefficients of the null model
end

%% Compute p-values
for i = 1:max_degree
    n = null_rich_coefs(:, i);   % Nulls at that degree threshold
    e = consensus_rich_coefs(i); % Empirical at that degree threshold
    
    % Compute how many null model rich club coefficients are greater than the consensus network's rich club coefficient
    p(i) = sum(n >= e) / null_sample;
end

%% Calculate normalized rich club coefficients
mean_rich_coefs = mean(null_rich_coefs);                   % Mean rich club coefficients from the null models
norm_rich_coefs = consensus_rich_coefs ./ mean_rich_coefs; % Take observed rich club coefficients and divide it by the average null model rich club coefficients

%% Create three plots to visualize output (Supplementary Materials Figure 8B-D)

degree = degrees_und(consensus_network);  % Calculate observed degree for histogram

figure(1); clf(1);
sgtitle('Consensus Network');

% Subplot 1: Plot normalized rich club coefficients with '1' cut off line
subplot(1, 3, 1);
hold on
yline(1);
plot(norm_rich_coefs);
hold off
ylabel('Normalized Rich Club Coef');
xlabel('Degree Threshold');
b = gca;
b.FontSize = 15;
b.TickDir = 'out';
b.FontName = 'Arial';
box off;

% Subplot 2: Plot the p-values across degree thresholds, add red dot where p-value is significant
subplot(1, 3, 2);
hold on
yline(0.05, '--');
plot(p);
scatter(find(p < 0.05), p(find(p < 0.05)), 15, 'filled', 'r');
hold off
ylabel('P-value');
xlabel('Degree Threshold');
b = gca;
b.FontSize = 15;
b.TickDir = 'out';
b.FontName = 'Arial';
box off;

% Subplot 3: Plot a histogram of observed degree for the consensus network, make all degrees above the rich club cut-off red
subplot(1, 3, 3);
hold on
histogram(degree);
histogram(degree(find(degree >= 12)), 'FaceColor', 'red');  % Change '12' to the degree threshold determined by the p-values
hold off
ylabel('Frequency');
xlabel('Degree');
b = gca;
b.FontSize = 15;
b.TickDir = 'out';
b.FontName = 'Arial';
box off;

%% Define rich club nodes

% Find nodes with the a degree at and above the cut-off (change to correct threshold)
rich_club_nodes = find(degree>=12==1); 

