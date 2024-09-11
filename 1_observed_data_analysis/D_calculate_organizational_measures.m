%% Calculate organizational measures
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script takes binarized networks and calculates local and global
organizational measures using the Brain Connectivity Toolbox
(https://sites.google.com/site/bctnet/home?authuser=0)

Toolbox publication:
Rubinov, M., & Sporns, O. (2010). Complex network measures of brain 
connectivity: uses and interpretations. Neuroimage, 52(3), 1059-1069.

%}

%% Add paths and load data
clear;clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Load binarized networks
load('binarized_connectomes.mat'); % or load('density_controlled_binarized_connectomes.mat') for density-controlled analysis
% Load atlas euclidean distance
load('euclidean_distance.mat');
% Load rich club nodes
load('rich_club_nodes.mat'); % or load('rich_club_nodes10.mat') for density-controlled analysis

nsub = size(binarized_connectomes,1); % Define number of participants


%% Calculate organizational measures

% Define names of chosen local and global measures as strings
localmeasures = string({'Degree','Betweenness','Clustering coefficient','Edge length','Local efficiency','Matching'});
globalmeasures = string({'Modularity','Characteristic path length','Global efficiency','Rich connections','Feeder connections','Local_connections',...
    'Rich connection lengths','Feeder connection lengths','Local connection lengths'});

% Initialize arrays
observed_local_statistics = zeros(nsub, 90, 6); 
observed_global_statistics = zeros(nsub, 9);   
            
% Loop over participants
for sub = 1:nsub
    % Extract the observed network for the current participant
    observed_network = squeeze(binarized_connectomes(sub,:,:));
    
    %%% Local statistics %%%
    % Degree
    observed_local_statistics(sub, :, 1) = degrees_und(observed_network);
    % Betweenness
    observed_local_statistics(sub, :, 2) = betweenness_bin(observed_network);
    % Clustering coefficient
    observed_local_statistics(sub, :, 3) = clustering_coef_bu(observed_network);
    % Edge length
    j = triu(observed_network, 1) > 0;                                     % Indices of upper triangular part of observed_network
    for node = 1:90
        h = j(node, :);                                                    % Edges from the node
        observed_local_statistics(sub, node, 4) = sum(edistance(node, h)); % Sum of distances from the node to other nodes
    end
    % Local efficiency
    observed_local_statistics(sub, :, 5) = efficiency_bin(observed_network, 1);
    % Matching
    observed_local_statistics(sub, :, 6) = mean(matching_ind(observed_network) + matching_ind(observed_network)');

    %%% Global statistics %%%
    % Maximum modularity
    [~, observed_global_statistics(sub, 1)] = modularity_und(observed_network);
    % Characteristic path length and global efficiency
    dist = distance_bin(observed_network);
    [observed_global_statistics(sub, 2), observed_global_statistics(sub, 3)] = charpath(dist, 0, 0);
    
    % Define connection types and lengths
    clear c                          % Clear variable from preivous loop
    half = triu(observed_network);   % Only take half on the matrix
    [c(:,1) c(:,2)] = find(half==1); % Find connections (save index of connected nodes)

    % Initialize
    counts = zeros(length(c),3); 
    con_length = zeros(length(c),3); 
    
    % Loop through number of connections
    for i = 1:length(c)
        if ismember(c(i,1),rich_club_nodes) && ismember(c(i,2),rich_club_nodes)       %% If both nodes are rich nodes = rich
            counts(i,1) = 1;
        elseif ismember(c(i,1),rich_club_nodes) || ismember(c(i,2),rich_club_nodes)   %% If one node is a rich club = feeder
            counts(i,2) = 1;
        elseif ~ismember(c(i,1),rich_club_nodes) && ~ismember(c(i,2),rich_club_nodes) %% If neither node is a rich club = local
            counts(i,3) = 1;
        else
            disp(sprintf('Missing condition for connention between nodes %g and %g...',c(i,1),c(i,2)));
        end
        % Find length of the connection 
        con_length(i,find(counts(i,:)==1)) = edistance(c(i,1),c(i,2)); 
    end
    % Save total number of connection types
    observed_global_statistics(sub,4) = nnz(counts(:,1)); % Rich connections
    observed_global_statistics(sub,5) = nnz(counts(:,2)); % Feeder connections
    observed_global_statistics(sub,6) = nnz(counts(:,3)); % Local connections
    % Save average length of connection types
    observed_global_statistics(sub,7) = mean(con_length(find(con_length(:,1)>0),1));
    observed_global_statistics(sub,8) = mean(con_length(find(con_length(:,2)>0),2));
    observed_global_statistics(sub,9) = mean(con_length(find(con_length(:,3)>0),3));
end
