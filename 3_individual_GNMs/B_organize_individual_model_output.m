%% Organize individual generative model outputs
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script takes individual generative model output files and puts them
into one struct to be used in further analysis.

%}

%% Add paths and load data
clear;clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Set save directory 
sdir = '/set/your/path/';               % <<<<<<<<<< SET

% Load binarized networks
load('binarized_connectomes.mat');
% Set number of participants and number of generative runs
nsub = size(binarized_connectomes,1); 
nruns = 25;                        

%% Load generative model data and put it into one struct

% Initialize
energy_sample = zeros(nsub,nruns);
ks_sample = zeros(nsub,nruns,4);
networks_sample = cell(nsub,1);
parameters_sample = zeros(nsub,nruns,2);
errors = zeros(nsub,1);

% Loop over participants
for sub = 1:nsub
    try % Load this paticipant's generative models 
        load(sprintf('generative_output_sub-%g.mat',sub));
    catch
        % Keep if it doesn't load
        errors(sub) = 1;
        % Print missing participant
        disp(sprintf('Participant %g non-existant',sub));
    end
    % Save variables
    energy_sample(sub,:,:) = output.energy;
    ks_sample(sub,:,:) = output.ks;
    networks_sample{sub} = output.networks;
    parameters_sample(sub,:,:) = output.params;
    
    % Clear the variable
    clear output
    
    % Display
    disp(sprintf('Participant %g loaded',sub));
end

% Save data as a struct
individual_generative_models = struct;
individual_generative_models.energy = energy_sample;
individual_generative_models.ks = ks_sample;
individual_generative_models.networks = networks_sample;
individual_generative_models.parameters = parameters_sample;
individual_generative_models.procedure = string({'Grid search n=10000 parameters eta [-3 0] and gamma [0.1 0.6] limits'}); % Include descriptives of how the models were run
save('individual_generative_models.mat','individual_generative_models','-v7.3');