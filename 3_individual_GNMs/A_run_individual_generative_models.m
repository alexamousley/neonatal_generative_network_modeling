%% Run individual generative models 
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script takes individual networks and runs the matching generative
models.

%}

%% Add paths and load data
clear;clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Set save directory 
sdir = '/set/your/path/';               % <<<<<<<<<< SET

% Load binarized networks
load('binarized_connectomes.mat');
% Load atlas euclidean distances
load('euclidean_distance.mat');

nsub = size(binarized_connectomes,1); % Define number of participants

%% Initialise model information
% Define model type
modeltype = string({'matching'});

% Set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'}, {'powerlaw'}];

%% Run grid search for parameter selection
% Choose upper and lower limits of parameters
eta = [-3, 0];    % Eta = 'cost' parameter
gam = [0.1, 0.6]; % Gamma = 'value' parameter

% Choose how many simulations to run (= how many parameter pairs)
nruns = 2;    

% Run grid search
[p, q] = meshgrid(linspace(eta(1), eta(2), sqrt(nruns)),...
    linspace(gam(1), gam(2), sqrt(nruns)));
params = unique([p(:) q(:)], 'rows');  % Create a list of parameter pairings
nparams = size(params, 1);             % Number of parameter pairings (equal to nruns)


%% Run generative models

% Loop through every participant
for sub = 1:nsub
    % select target connectome
    observed_network = squeeze(binarized_connectomes(sub,:,:)); 

    % Initialize X and calculate organizational measures of the observed network 
    x = cell(4,1);
    x{1} = sum(observed_network,2);
    x{2} = clustering_coef_bu(observed_network);
    x{3} = betweenness_bin(observed_network)';
    x{4} = edistance(triu(observed_network,1) > 0);

    % Calculate seed network 
    proportion    = 0.95;                                            % Set proportion of edges in common across participants
    connections   = squeeze(mean(binarized_connectomes,1));  % Define connections from the average network
    index         = find(connections>proportion);                    % Find connections common across the set proportion of participants
    A             = zeros(size(connections));                        % Create empty array for seed network
    A(index)      = 1; A = upper(A);                                 % Create seed network (add 1 at every index)

    % Define the number of connections in the observed network (the generative models will stop when they reach the same number of connections)
    m    = nnz(observed_network)/2;  
    % Define number of nodes in the network
    n    = length(observed_network);  
    
    % Print text
    disp(sprintf('Running participant %g using %s model...',sub,modeltype));
    
    % Initialize
    output = struct;
    output.params = zeros(nparams,2);
    output.energy = zeros(nparams);
    output.ks = zeros(nparams,4);
    output.networks = zeros(m,nparams);  
    
    % Run model
    B = generative_model(A,edistance,m,modeltype,modelvar,params);
    
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
    output.params = params;
    output.energy = E;
    output.ks = K;
    output.networks = B;
    
    % Change directory
    cd(sdir);
    % Save pariticpant's generative models
    save(sprintf('generative_output_sub-%g.mat',sub),'output','-v7.3');
end

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