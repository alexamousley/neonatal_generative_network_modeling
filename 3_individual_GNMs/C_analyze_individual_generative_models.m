%% Run individual generative models 
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script explores model fits and orgnaization of individual generative
models.

(1) Identify the best fit models based on lowest energy
(2) Create mean energy landscape
(3) Calculate organizational measures for the simulated networks
(4) Explore spatial embedding via the cumulative probability disribution
(5) Calculate Topological fingerprints and dissimilarity

%}

%% Add paths and load data
clear;clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Set directory to save data
datadir = '/set/your/path';             % <<<<<<<<<< SET

% Load binarized networks
load('binarized_connectomes.mat');
% Load individual generative models stuct
load('individual_generative_models.mat');
% Load atlas euclidean distances
load('euclidean_distance.mat');
% Load birth age
load('GA_at_birth.mat');
term_index = GA_at_birth >= 37; % Define which participants were born at term
% Load observed organizational measures
load('observed_global_statistics.mat')
load('observed_local_statistics.mat');
% Load rich club nodes
load('rich_club_nodes.mat');

% Set number of participants and generative models
nsub = size(binarized_connectomes,1); 
nruns = length(individual_generative_models.energy);
%% (1) Identify best fitting parameters and lowest energy model

% Initialize
params_sorted = zeros(nsub,nruns,2);
Asynth_sorted = cell(size(individual_generative_models.networks)); 
Eall_sorted = zeros(size(individual_generative_models.energy));

% Loop over each participant
for sub = 1:nsub
    h = squeeze(individual_generative_models.energy(sub,:));                    % Select energy for subject for model type
    [~,i] = sort(h);                                                 % Create index of acending order of energy
    Eall_sorted(sub,:) = individual_generative_models.energy(sub,i);            % Sort energy
    params_sorted(sub,:,:) = individual_generative_models.parameters(sub,i,:);  % Sort parameters by energy index
    Asynth_sorted{sub} = individual_generative_models.networks{sub}(:,i);       % Sort generative models by energy index
    disp(sprintf('For subject %g, the lowest energy is %g (achieved with parameters p = %g and q = %g)',sub,round(Eall_sorted(sub,1),2),params_sorted(sub,1,1),params_sorted(sub,1,2)));
end

% Identify the best-fitting parameters by taking the average of the top 10 networks
top = 10;
for sub = 1:nsub
    % Calculate means
    mean_eta = mean(squeeze(params_sorted(sub,1:top,1)));
    mean_gamma = mean(squeeze(params_sorted(sub,1:top,2)));
    mean_top_params(sub,:) = [mean_eta mean_gamma];
end

%% (2) Plot energy landscapes

% Set up variables
emean = mean(individual_generative_models.energy);           % Take mean energy 
p = squeeze(individual_generative_models.parameters(1,:,:)); % Define parameters

% Plot energy landscape (Figure 4C)
h=figure(2); h.Position=[100 100 600 400];
hold on
scatter(p(:,1),p(:,2),2000,emean,'.');                                     % Plot average energies
s = scatter(mean_top_params(:,1),mean_top_params(:,2),20,'white','filled');% Plot parameter location of individual's best-fit models
hold off
colormap(gca,'hot');
xlim([-3 0]);
ylim([0.1 0.6]);
xlabel('Eta'); ylabel('Gamma');
set(gca,'FontSize',40);
b = gca; b.TickDir = 'out'; b.FontName='Arial'; box off;


%% (3) Calculate organizational measures for simulated networks

% Define names of chosen local and global measures as strings
localmeasures = string({'Degree','Betweenness','Clustering coefficient','Edge length','Local efficiency','Matching'});
globalmeasures = string({'Modularity','Characteristic path length','Global efficiency','Rich connections','Feeder connections','Local_connections',...
    'Rich connection lengths','Feeder connection lengths','Local connection lengths'});

% Initialize
simulated_local_statistics  = zeros(nsub,90,6);
simulated_global_statistics = zeros(nsub,9);

% Loop over participants
for sub  = 1:nsub

    % Take simulated network
    s = squeeze(Asynth_sorted{sub});
    % Reformat networks into a region x region matrix
    A = zeros(90);
    A(s(:,1)) = 1;
    simulated_network = A + A';  

    %%% Local statistics %%%
    % Degree
    simulated_local_statistics(sub,:,1) = degrees_und(simulated_network);
    % Betweeness
    simulated_local_statistics(sub,:,2) = betweenness_bin(simulated_network);
    % Clustering
    simulated_local_statistics(sub,:,3) = clustering_coef_bu(simulated_network);
    % Edge length 
    r = triu(simulated_network,1)>0;        % Indices of upper triangular of network
    for node = 1:90                 
        s = r(node,:);                      % Edges from node   
        simulated_local_statistics(sub,node,4) = ...
            sum(edistance(node,s));         % Sum of distances from the node to other nodes
    end
    % Local efficiency
    simulated_local_statistics(sub,:,5) = efficiency_bin(simulated_network,1);
    % Matching
    simulated_local_statistics(sub,:,6) = mean(matching_ind(simulated_network)+matching_ind(simulated_network)');
    
    %%% Global statistics %%%
    % Modularity 
    [~,simulated_global_statistics(sub,1)] = modularity_und(simulated_network);
    % Characteristic path length and global efficiency
    dist = distance_bin(simulated_network);
    [simulated_global_statistics(sub,2) simulated_global_statistics(sub,3)] = charpath(dist,0,0);
    
    % Define connection types and lengths
    clear c                          % Clear variable from preivous loop
    half = triu(simulated_network);   % Only take half on the matrix
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
    simulated_global_statistics(sub,4) = nnz(counts(:,1)); % Rich connections
    simulated_global_statistics(sub,5) = nnz(counts(:,2)); % Feeder connections
    simulated_global_statistics(sub,6) = nnz(counts(:,3)); % Local connections
    % Save average length of connection types
    simulated_global_statistics(sub,7) = mean(con_length(find(con_length(:,1)>0),1));
    simulated_global_statistics(sub,8) = mean(con_length(find(con_length(:,2)>0),2));
    simulated_global_statistics(sub,9) = mean(con_length(find(con_length(:,3)>0),3));
end
%% (4) Visualize cumulative probability disribution (Figure 4F) 
% CPD plots were adapted from Danyal Akarca's code (danyal.akarca@mrc-cbu.cam.ac.uk)

% Choose which measure to plot
measure = 1;

% Take the statistics
obs = mean(squeeze(observed_local_statistics(:,:,measure)))';
sim = mean(squeeze(simulated_local_statistics(:,:,measure)))';

% Make common bins
binEdges =  [-inf;sort([obs;sim]);inf];

% Compute the bin counts
binCounts1 = histc(obs,binEdges,1);
binCounts2 = histc(sim,binEdges,1);

% Compute a cumulative probability
sumCounts1 = cumsum(binCounts1)./sum(binCounts1);
sumCounts2 = cumsum(binCounts2)./sum(binCounts2);

% Compute the cumulative density
sampleCDF1 = sumCounts1(1:end-1);
sampleCDF2 = sumCounts2(1:end-1);

% Plot (Figure 4F and Supplementary Figure 6)
figure(1); clf;
plot(smooth(sampleCDF1),'linewidth',6,'color',[0 0 0]);
hold on;
plot(smooth(sampleCDF2),'linewidth',6,'color',[0.5977 0.1953 0.7969]);
set(gca,'TickLength',[0 0]);
title(''); xticks([0 200]); yticks([0 1]);
set(gca,'color','none'); 
xlabel(localmeasures(measure)); ylabel(sprintf('F(%s)',localmeasures(measure)));
l=legend(string({'Observed','Simulated'}),'location','northeastoutside');
set(l,'color','none');
set(gca,'FontSize',20);
b = gca; b.TickDir = 'out'; b.FontName='Arial'; box off;

%% (5) Topological fingerprints and dissimilarity

% Initialize
observed_tf = zeros(nsub,6,6);
simulated_tf = zeros(nsub,6,6);
tfdissimilarity = zeros(nsub,6);

% Loop through participants
for sub = 1:nsub
    % Take participant's statistics 
    observed = squeeze(observed_local_statistics(sub,:,:));
    simulated = squeeze(simulated_local_statistics(sub,:,:));
    
    % Remove non-cortical and all 0 nodes
    observed(71:80,:) = [];
    simulated(71:80,:) = [];
    
    % Calculate topological fingerprints
    observed_tf(sub,:,:) = corr(observed);
    simulated_tf(sub,:,:) = corr(simulated);
    
    % Calculate topoloical dissimilarity
    tfdissimilarity(sub,:) = norm(squeeze(observed_tf(sub,:,:))-squeeze(simulated_tf(sub,:,:)));
end


% Plot TF matrix (Figure 4C)

% Define average topological fingerprints
mean_observed_tf = squeeze(mean(observed_tf)); 
mean_simulated_tf = squeeze(mean(simulated_tf));

% Plot
figure(25);clf(25); set(gcf,'Position',[100 1300 1000 1000]);
subplot(1,2,1); 
imagesc(mean_observed_tf);
subtitle('Observed');
set(gca,'xtick',1:6,'xticklabel',localmeasures,'fontsize',15);
set(gca,'ytick',1:6,'yticklabel',localmeasures,'fontsize',20);
xtickangle(40); set(gca,'FontSize',20); 
b = gca; b.TickDir = 'out'; b.FontName='Arial'; box off;
subplot(1,2,2); 
imagesc(mean_simulated_tf);
subtitle('Simulated'); set(gca,'FontSize',20);
set(gca,'xtick',[],'xticklabel',[]);
set(gca,'ytick',[],'yticklabel',[]); colorbar;
b = gca; b.TickDir = 'out'; b.FontName='Arial'; box off;

