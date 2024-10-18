%% Analyze initial generative models 
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script explores the output from 13 different generative models fit to
the consensus network.

(1) Identify the best fit models based on lowest energy
(2) Create mean energy landscapes
(3) Calculate organizational measures for the simulated networks
(4) Explore spatial embedding of local organizational structure
(5) Calculate Topological fingerprints and dissimilarity

%}

%% Add path and load data 
clear;clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Load binarized networks
load('binarized_connectomes.mat');
% Load consensus network
load('consensus_network.mat');
% Load initial generative models
load('initial_generative_models.mat');
% Load atlas euclidean distance
load('euclidean_distance.mat');

%% Define initial variables

% Define model types (13 different 'values')
modeltype = string({'sptl', 'neighbors', 'matching', 'clu-avg', 'clu-min', 'clu-max', 'clu-diff', 'clu-prod', 'deg-avg', 'deg-min', 'deg-max', 'deg-diff', 'deg-prod'});
nmodels = length(modeltype);  % number of models
% Display full names of models, just for visualizations
modelname = string({'Spatial','Neighbors','Matching',...
    'Clustering Average','Clustering Minimum','Clusetering Maximum','Clustering Difference','Clustering Product',...
    'Degree Average','Degree Minimum','Degree Maximum','Degree Difference','Degree Product'});

% Set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'}, {'powerlaw'}];

% Define KS statistics
ks_names = ({'Degree','Clustering','Betweenness','Edge Length'});

% Define number of runs
nruns = length(initial_generative_models.params);
% Define size of network
nsize = height(initial_generative_models.networks{1}); 

%% (1) Identify best fitting models by finding the lowest energy and associated parameters

% Initialize arrays (the first element of 'sorted' arrays will be the best fit model information)
params_sorted = zeros(nmodels,nruns,2);     % Parameters
Asynth_sorted = zeros(nmodels,nsize,nruns); % Generative Models
Eall_sorted = zeros(size(initial_generative_models.energy));   % Energies

for model = 1:nmodels
    modelE = squeeze(initial_generative_models.energy(model,:));                    % Select energy for model type
    [~,modelI] = sort(modelE);                                                      % Create index of that sorts energy (from lowest to highest)
    Eall_sorted(model,:) = initial_generative_models.energy(model,modelI);          % Create a new variable of energy that is sorted 
    params_sorted(model,:,:) = initial_generative_models.params(model,modelI,:,:);  % Sort parameters by energy index
    model_networks = squeeze(initial_generative_models.networks{model});            % Pull which model to sort 
    Asynth_sorted(model,:,:) = model_networks(:,modelI);                            % Sort generative models by energy index
end 

% Print descriptives for each model type 
% Initialize
elow_per_model = [];  % Energy
eta_per_model = [];   % Eta
gamma_per_model = []; % Gamma

% Loop through models
for model = 1:nmodels
    elow_per_model(model,:) = Eall_sorted(model,1);      % Lowest energy
    eta_per_model(model,:) = params_sorted(model,1,1);   % Best eta
    gamma_per_model(model,:) = params_sorted(model,1,2); % Best gamma
    % Energy mean and standard deviation
    eMean(model,:) = mean(Eall_sorted(model,:));
    eSD(model,:) = std(Eall_sorted(model,:));
end  

% Display means and standard deviations for each model
networktable = table(modelname',eMean,elow_per_model,eSD,eta_per_model,gamma_per_model,'VariableNames',...
    {'Rule','Mean Energy','Lowest Energy','Std Energy','Best Eta','Best Gamma '});
disp(networktable);

%% (2) Plot energy landscape for each model type

% Plot all model's energy landscapes (Supplementary materials Figure 5B)
figure(1);clf(1);
for model = 1:nmodels
    % Take the model's energy
    e = squeeze(initial_generative_models.energy(model,:));
    % Create landscape plot
    subplot(4,4,model);
    scatter(initial_generative_models.params(model,:,1),initial_generative_models.params(model,:,2),400,e,'.');
    t = title(modelname(model));
    t.FontSize = 10;
    xlabel('Eta'); ylabel('Gamma'); 
    caxis([0 1]); c = colorbar; c.Label.String = 'Energy';
    b = gca; b.TickDir = 'out';
    colormap(gca,'hot');
end

% Plot landscape for one model type (Figure 4A)
model = 3; % Choose model
figure(2);clf(2);
hold on
e = squeeze(initial_generative_models.energy(model,:));
scatter(initial_generative_models.params(model,:,1),initial_generative_models.params(model,:,2),100,e,'.'); 
scatter(params_sorted(model,1,1),params_sorted(model,1,2),'g','filled');
hold off
xticks([-8 0]); yticks([-8 8]); title(modelname(model)); xlabel('Eta'); ylabel('Gamma'); 
caxis([0 1]); c = colorbar; c.Label.String = 'Energy'; colormap(gca,'hot');
b = gca; b.TickDir = 'out';
set(gca, 'FontName', 'Arial');       
c.Label.FontName = 'Arial';  

%% (3) Calculate organizational measures for observed and simulated networks

% Define names of chosen local and global measures as strings
localmeasures = string({'Degree','Betweenness','Clustering coefficient','Edge length','Local efficiency','Matching'});
globalmeasures = string({'Modularity','Characteristic path length','Global efficiency'});

% Initialize
observed_local_statistics   = zeros(90,6);
simulated_local_statistics  = zeros(nmodels,90,6);
observed_global_statistics  = zeros(1,3);
simulated_global_statistics = zeros(nmodels,3);

% Loop over models
for model = 1:nmodels
    
    % Take observed network (which is the consensus network)
    observed_network = consensus_network; 
    
    % Take simulated network
    s = squeeze(Asynth_sorted(model,:,:));
    % Reformat networks into a region x region matrix
    A = zeros(90);
    A(s(:,1)) = 1;
    simulated_network = A + A';  

    %%% Local statistics %%%
    % Degree
    observed_local_statistics(:,1) = degrees_und(observed_network);
    simulated_local_statistics(model,:,1) = degrees_und(simulated_network);
    % Betweeness
    observed_local_statistics(:,2) = betweenness_bin(observed_network);
    simulated_local_statistics(model,:,2) = betweenness_bin(simulated_network);
    % Clustering
    observed_local_statistics(:,3) = clustering_coef_bu(observed_network);
    simulated_local_statistics(model,:,3) = clustering_coef_bu(simulated_network);
    % Edge length
    j = triu(observed_network,1)>0;     % Indices of upper triangular of network
    r = triu(simulated_network,1)>0;
    for node = 1:90
        h = j(node,:);                  % Edges from node
        s = r(node,:);
        observed_local_statistics(node,4)  = ...
            sum(edistance(node,h));     % Sum of distances from the node to other nodes
        simulated_local_statistics(model,node,4) = ...
            sum(edistance(node,s));
    end
    % Local efficiency
    observed_local_statistics(:,5) = efficiency_bin(observed_network,1);
    simulated_local_statistics(model,:,5) = efficiency_bin(simulated_network,1);
    % Matching
    observed_local_statistics(:,6) = mean(matching_ind(observed_network)+matching_ind(observed_network)');
    simulated_local_statistics(model,:,6) = mean(matching_ind(simulated_network)+matching_ind(simulated_network)');
    
    %%% Global statistics %%%
    % Modularity 
    [~,observed_global_statistics(1)] = modularity_und(observed_network);
    [~,simulated_global_statistics(model,1)] = modularity_und(simulated_network);
    % Characteristic path length and global efficiency
    dist = distance_bin(observed_network);
    [observed_global_statistics(2) observed_global_statistics(3)] = charpath(dist,0,0);
    dist = distance_bin(simulated_network);
    [simulated_global_statistics(model,2) simulated_global_statistics(model,3)] = charpath(dist,0,0);
end

%% (4) Calculate spatial embedding (the correlation between observed and simulated network's local organization)

% Initialize
spatial_embed = zeros(nmodels,length(modeltype),2); % Save r as first element and p-value as second
% Loop over models
for model = 1:nmodels
    % Loop over number of local statistics
    for measure = 1:length(localmeasures)
        % Take measure for observed and simulated networks
        observed = squeeze(observed_local_statistics(:,measure));
        simulated = squeeze(simulated_local_statistics(model,:,measure));
        
        % Remove non-cortical and all 0 nodes
        observed(71:80,:) = [];
        simulated(:,71:80) = [];
        
        % Calcuate correlation
        [r p] = corr(observed,simulated');
        spatial_embed(model,measure,:) = [r p];
    end
end

% Plot spatial embedding for one model type (Supplementary materials Figure 5D)

% Choose model type
model = 3;
for measure = 1:6
    % Take networks
    observed = squeeze(observed_local_statistics(:,measure));
    simulated = squeeze(simulated_local_statistics(model,:,measure));
    
    % Define correlation values
    r = spatial_embed(model,measure,1);
    p = spatial_embed(model,measure,2);

    % Compute line of best fit
    a = polyfit(observed,simulated,1); 
    f = polyval(a,observed);
    
    % Plot
    figure(3);
    subplot(2,3,measure);
    scatter(observed,simulated,'filled','k','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2);
    hold on;
    a = plot(observed,f,'k'); a.LineWidth = 5;
    % Change color of correlation line if p is significant and print values as subtitle
    if round(p,4) < 0.001 
        a.Color = 'red';
        s = subtitle(sprintf('r = %g, p < 0.001',round(r,3))); s.FontSize =13;
    end 
    if round(p,4) >= 0.001 && round(p,4) < 0.05 
        a.Color = 'blue';
        s = subtitle(sprintf('r = %g, p < 0.05',round(r,3))); s.FontSize =13;
    end 
    
    t = title(sprintf('%s',localmeasures(measure))); t.FontSize = 15;
    ylim([min(simulated),max(simulated)]);
    xlabel('Observed','FontSize',20); ylabel('Simulated','FontSize',20);
end

%% (5) Calculate topological fingerprints (the correlation between local organizational measures)

% Initialize
observed_tf = zeros(6,6);           % Correlation marix for consensus network
simulated_tf = zeros(nmodels,6,6);  % Correlation matrices for each model
tfdissimilarity = zeros(nmodels,6); % Topological dissmilarity

% Loop through each model
for model = 1:nmodels
    % Take observed and simulated statistics
    observed = observed_local_statistics;
    simulated = squeeze(simulated_local_statistics(model,:,:));
    
    % Remove non-cortical and all 0 nodes
    observed(71:80,:) = [];
    simulated(71:80,:) = [];
    
    % Calculate topological fingerprints
    observed_tf = corr(observed);
    simulated_tf(model,:,:) = corr(simulated);
    
    % Calculate topoloical dissimilarity
    tfdissimilarity(model,:) = norm(observed_tf-squeeze(simulated_tf(model,:,:)));
end

% Plot TF matrix (Supplementary Figure 5A)
matrix = observed_tf; % Choose which topological fingerprint to plot

% Plot
figure(25);clf(25); set(gcf,'Position',[100 1000 900 900]);
imagesc(matrix);
set(gca,'xtick',[1:6],'xticklabel',localmeasures,'fontsize',15);
set(gca,'ytick',[1:6],'yticklabel',localmeasures,'fontsize',20);
xtickangle(40); set(gca,'FontSize',20); colorbar;
b = gca; b.TickDir = 'out'; b.FontName='Arial'; box off;

