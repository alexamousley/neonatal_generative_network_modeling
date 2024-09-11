%% Run initial generative models 
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script takes the consensus network and runs 13 generative 
models with different 'value' metrics. The generative models are created 
from the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/home?authuser=0)
function. The parameters are selected using a grid search and the seed network 
is defined by connections in common across a pre-defined percent (e.g., 95%) 
of the whole sample's binarized networks.

Note: These generative models are probabilistic, which means that if you run 
your own simulations, they will not be identical to the simulations we have 
published.

%}

%% Add paths and load data
clear;clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Set save directory 
sdir = '/set/your/path/';               % <<<<<<<<<< SET

% Load binarized networks
load('binarized_connectomes.mat');
% Load consensus network
load('consensus_network.mat');
% Load atlas euclidean distance
load('euclidean_distance.mat');

nsub = size(binarized_connectomes,1); % Set number of participants      

%% Initialise model information
% Define model types (13 different 'values')
modeltype = string({'sptl', 'neighbors', 'matching', 'clu-avg', 'clu-min', 'clu-max', 'clu-diff', 'clu-prod', 'deg-avg', 'deg-min', 'deg-max', 'deg-diff', 'deg-prod'});
nmodels = length(modeltype);  % number of models

% Set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'}, {'powerlaw'}];

%% Run grid search for parameter selection
% Choose upper and lower limits of parameters
eta = [-8, 0]; % Eta = 'cost' parameter
gam = [-8, 8]; % Gamma = 'value' parameter

% Choose how many simulations to run (= how many parameter pairs)
nruns = 25;    

% Run grid search
[p, q] = meshgrid(linspace(eta(1), eta(2), sqrt(nruns)),...
    linspace(gam(1), gam(2), sqrt(nruns)));

params = unique([p(:) q(:)], 'rows');  % Create a list of parameter pairings
nparams = size(params, 1);             % Number of parameter pairings (equal to nruns)

     
%% Run generative models

% Initialize X and calculate organizational measures of the consensus network for model-fitting (energy equation)
x = cell(4,1);
x{1} = sum(consensus_network,2);
x{2} = clustering_coef_bu(consensus_network);
x{3} = betweenness_bin(consensus_network)';
x{4} = edistance(triu(consensus_network,1) > 0);

% Calculate seed network 
proportion    = 0.95;                                            % Set proportion of edges in common across participants
connections   = squeeze(mean(binarized_connectomes,1));  % Define connections from the average network
index         = find(connections>proportion);                    % Find connections common across the set proportion of participants
A             = zeros(size(connections));                        % Create empty array for seed network
A(index)      = 1; A = upper(A);                                 % Create seed network (add 1 at every index)

% Define the number of connections in the consensus network (the generative models will stop when they reach the same number of connections)
m    = nnz(consensus_network)/2;  
% Define number of nodes in the network
n    = length(consensus_network);      

% Initialize struct to store generative models
initial_generative_models = struct;
initial_generative_models.params = zeros(nmodels,nparams,2);
initial_generative_models.energy = zeros(nmodels,nparams);
initial_generative_models.ks = zeros(nmodels,nparams,4);
initial_generative_models.networks = zeros(nmodels,m,nparams);   

% Loop through different model types
for model = 1:nmodels
    % Print text
    disp(sprintf('Running consensus network with %s model...',modeltype(model)));
    
    % Run the model
    B = generative_model(A,edistance,m,modeltype{model},modelvar,params);
    
    % Initialize array for KS statistics
    K = zeros(nparams,4);
    % Loop through runs (different parameters)
    for iB = 1:nparams   
        
        % Turn the generative model B into a region x region matrix
        b = zeros(n);
        b(B(:,iB)) = 1;
        b = b + b';  
        
        % Initialize Y and calculate organizational measures of the generative model for model-fitting (energy equation)
        y = cell(4,1);
        y{1} = sum(b,2);
        y{2} = clustering_coef_bu(b);
        y{3} = betweenness_bin(b)';
        y{4} = edistance(triu(b,1) > 0);
        
        % Calculate the KS statistics
        for j = 1:4
            K(iB,j) = fcn_ks(x{j},y{j});
        end
    end
    
    % Define the output of the energy equation (the maximum KS statistic)
    E = max(K,[],2);
    
    % Store the output
    initial_generative_models.params(model,:,:) = params;
    initial_generative_models.energy(model,:) = E;
    initial_generative_models.ks(model,:,:) = K;
    initial_generative_models.networks(model,:,:) = B;
end

%% Save generative model output
% Change directory
cd(sdir);
% Save file
save('initial_generative_models.mat','initial_generative_models','-v7.3'); % Likely will need the '-v7.3' due to size of the struct

%% Define KS function
% This function was written by Danyal Akarca (danyal.akarca@mrc-cbu.cam.ac.uk),
% as a part of his generative modelling script that can be found here:
% https://github.com/DanAkarca/generativenetworkmodel/blob/master/Scripts/iii.%20Running%20initial%20generative%20models.m

function kstat = fcn_ks(x1,x2)
    binEdges    =  [-inf ; sort([x1;x2]) ; inf];
    binCounts1  =  histc (x1 , binEdges, 1);
    binCounts2  =  histc (x2 , binEdges, 1);
    sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
    sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);
    sampleCDF1  =  sumCounts1(1:end-1);
    sampleCDF2  =  sumCounts2(1:end-1);
    deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
    kstat = max(deltaCDF);
end
