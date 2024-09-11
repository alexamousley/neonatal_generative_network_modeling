%% Create group developmental generative models 
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script uses two parameters, which are pre-defined to represent a
'term' and 'preterm' network, and runs 1000 iterations of the matching
generative model. The difference between these models and the previous
models is that we are saving the network at every iteration (i.e., every
time a new connection is added), in order to explore how the network
developes. This involved adapting the Brain Connectivity Toolbox's
generative_model.m function.


%}

%% Add path and load data 
clear;clc;


% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Set save directory 
sdir = '/set/your/path/';            % <<<<<<<<<< SET


% Load binarized networks  
load('binarized_connectomes.mat');
% Load consensus network
load('consensus_network.mat');
% Load initial generative models
load('initial_generative_models.mat');
% Load atlas euclidean distances
load('euclidean_distance.mat');

%% Set up variables

% Define parameters
params = [-1.7435, 0.3249;-1.8457, 0.3374]; % First pair is 'term' and second pair is 'preterm'
nparams = size(params,1);                 % Number of parameters

% Choose number of runs
nruns = 25;

% Define number of connections (based on consensus network)
m = nnz(consensus_network)/2;

% Set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'}, {'powerlaw'}];

% Minimum edge
epsilon = 1e-5;
%% Run generative models

% Initialize
developmental_gnms = struct;
developmental_gnms.params = params;                   % Save parameters
developmental_gnms.networks = zeros(2,nruns,m,90,90); % Save networks at every iteration (computationally expensive)

% Loop through number of runs
for run = 1:nruns
    % Loop through number of parameters
    for iparam = 1:nparams
        
        % Calculate seed network
        proportion    = 0.95;                                            % Set proportion of edges in common across participants
        connections   = squeeze(mean(binarized_connectomes,1));  % Define connections from the average network
        index         = find(connections>proportion);                    % Find connections common across the set proportion of participants
        A             = zeros(size(connections));                        % Create empty array for seed network
        A(index)      = 1; A = upper(A);                                 % Create seed network (add 1 at every index)

        % The code below is adapted from the Brain Connectivity Toolbox
        % 'generative_model.m function. The adaptation involved selecting
        % only the 'matching' model, as well as saving every iteration of
        % the model's development
        
        % Run matching model
        n           = length(edistance);
        b           = zeros(m,nparams);
        K           = matching_ind(A);
        K           = K + K';

        eta = params(iparam,1);
        gam = params(iparam,2);
        K = K + epsilon;
        n = length(edistance);
        mseed = nnz(A)/2;
        mv1 = modelvar{1};
        mv2 = modelvar{2};
        switch mv1
            case 'powerlaw'
                Fd = edistance.^eta;
            case 'exponential'
                Fd = exp(eta*edistance);
        end
        switch mv2
            case 'powerlaw'
                Fk = K.^gam;
            case 'exponential'
                Fk = exp(gam*K);
        end
        Ff = Fd.*Fk.*~A;
        [u,v] = find(triu(ones(n),1));
        indx = (v - 1)*n + u;
        P = Ff(indx);
        % Save the first parameterised K in Fk
        FKall(1,:,:)  = Fk;
        % Save the first probabilities
        Ff(isinf(Ff)) = 0;
        Pall(1,:,:)   = Ff;
        step = 2; 
        for ii = (mseed + 1):m
            C = [0; cumsum(P)];
            r = sum(rand*C(end) >= C);
            uu = u(r);
            vv = v(r);
            A(uu,vv) = 1;
            A(vv,uu) = 1;
            updateuu = find(A*A(:,uu));
            updateuu(updateuu == uu) = [];
            updateuu(updateuu == vv) = [];
            updatevv = find(A*A(:,vv));
            updatevv(updatevv == uu) = [];
            updatevv(updatevv == vv) = [];
            c1 = [A(:,uu)', A(uu,:)];
            for i = 1:length(updateuu)
                j = updateuu(i);
                c2 = [A(:,j)' A(j,:)];
                use = ~(~c1&~c2);
                use(uu) = 0;  use(uu+n) = 0;
                use(j) = 0;  use(j+n) = 0;
                ncon = sum(c1(use))+sum(c2(use));
                if (ncon==0)
                    K(uu,j) = epsilon;
                    K(j,uu) = epsilon;
                else
                    K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                    K(j,uu) = K(uu,j);
                end
            end
            c1 = [A(:,vv)', A(vv,:)];
            for i = 1:length(updatevv)
                j = updatevv(i);
                c2 = [A(:,j)' A(j,:)];
                use = ~(~c1&~c2);
                use(vv) = 0;  use(vv+n) = 0;
                use(j) = 0;  use(j+n) = 0;
                ncon = sum(c1(use))+sum(c2(use));
                if (ncon==0)
                    K(vv,j) = epsilon;
                    K(j,vv) = epsilon;
                else
                    K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                    K(j,vv) = K(vv,j);
                end
            end
            switch mv2
               case 'powerlaw'
                    Fk = K.^gam;
                case 'exponential'
                    Fk = exp(gam*K);
            end
            % Save network at this iteration
            developmental_gnms.networks(iparam,run,ii,:,:) = A;
            Ff = Fd.*Fk.*~A;
            P = Ff(indx);
            Ff(isinf(Ff))   = 0; 
            % Change the step
            step = step+1;  
        end
        b(:,iparam)         = find(triu(A,1));
        disp(sprintf('Model %g, run %g is finished',iparam,run))
    end
end
 
%% Save models

% Change directory
cd(sdir);
% Save file
save('developmental_generative_models.mat','developmental_gnms','-v7.3'); % Likely will need the '-v7.3' due to size of the struct
