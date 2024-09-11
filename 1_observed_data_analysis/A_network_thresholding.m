%% Network Thresholding
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

This script takes raw networks, performs consensus thresholding followed by 
absolute thresholding and binarization. With thresholds as is (consesus =
0.6 and binarisation = 325), the resulting 'binarized_connectomes' variable
should be exactly equal to that of the 'binarized_connectomes' variable 
published.

Note: If you are using the example data, the consensus thresholding will 
need to be much lower that 60% as the data is random and therefore not many 
'connections' will be in common across the fake participants

%}
%% Add paths and load data
clear;clc;

% Add paths
run('/set/your/path/.../set_paths.m');  % <<<<< Add path to set_paths file

% Load unthresholded networks
load('unthresholded_connectomes.mat'); 

nsub = size(unthresholded_connectomes,1); % Set number of participants

%% Perform consensus thresholding
% Set threshold
set = 0.6;                     % Percentage threshold
threshold = floor(nsub * set); % Calculate the threshold value based on the number of participants

k = unthresholded_connectomes ~= 0;     % Find nonzero elements (essentially creating a mask)
u = squeeze(sum(k, 1));        % Sum the values along the first dimension, squeezing the result to remove singleton dimensions

% Keep/Remove indices
ind = u < threshold;           % Find indices where the sum is less than the threshold
indkeep = u > threshold;       % Find indices where the sum is greater than the threshold

% Apply Threshold
for sub = 1:nsub
    A = squeeze(unthresholded_connectomes(sub,:,:)); % Select one participant
    A(ind)=0;                               % Remove edges with index
    consensus_thresholded(sub,:,:)=A;       % Save network 
end

% Look at the number of connections before and after thresholding
for sub = 1:nsub
    % Get the original network (A) and thresholded network (B)
    A = squeeze(unthresholded_connectomes(sub,:,:));
    B = squeeze(consensus_thresholded(sub,:,:));
    
    % Calculate the number of connections before and after thresholding
    before_threshold(sub) = nnz(A)/2;    % Count the number of nonzero elements in A and divide by 2 (due to symmetric connections)
    after_threshold(sub) = nnz(B)/2;     % Repeat for B
end

%% Perform absolute thresholding and binarize networks

% Set threshold
thr = 325;

% Initialize variables
binarized_connectomes = [];    % Empty array to store binarised connectomes
density = [];                  % Empty array to store density values

% Apply thresholding
for sub = 1:nsub
    W = squeeze(consensus_thresholded(sub,:,:));       % Extract the consensus thresholded connectome for the current participant
    weighted_connectome = threshold_absolute(W, thr);  % Apply the absolute threshold to the connectome
    density(sub) = density_und(weighted_connectome);   % Calculate the density of the weighted network
    B = zeros(size(weighted_connectome));              % Create a matrix of zeros with the same size as the weighted network
    B(find(weighted_connectome)) = 1;                  % Find the connections above the threshold and set them to 1
    binarized_connectomes(sub,:,:) = B;                % Store the binarized connectome 
end

% Print results
disp(sprintf('At this threshold, a mean density of %g%% is produced across the sample.',100*mean(density)));

%% Perform density-controlled thresholding for the density-controlled analysis

thr  = [1:4000];     % Create various thresholds to iterate through
nthr = length(thr);  % Number of thresholds to try

% Initialize variables
binarized_connectomes = [];    % Empty array to store binarised connectomes
density = [];                  % Empty array to store density values

% Iterate through participants
for sub = 1:nsub 
    % Loop through each threshold option
    for t = 1:nthr 
        over = 0;                                            % Create variable for if density is over-shot (<10%)
        W     = squeeze(unthresholded_connectomes(sub,:,:)); % Extract the connectome for the current participant
        weighted_connectome = threshold_absolute(W, thr(t)); % Apply the absolute threshold to the connectome
        d = density_und(weighted_connectome);                % Calculate density
        if round(d,3) == 0.100 || over == 1                  % If the density is equal to 10% or density has been over-shot (indicated by over == 1)
            B = zeros(size(weighted_connectome));            % Create a matrix of zeros with the same size as the weighted network
            B(find(weighted_connectome)) = 1;                % Find the connections above the threshold and set them to 1
            binarized_connectomes(sub,:,:) = B;              % Store the binarized connectome 
            density(sub) = d;                                % Save density
            %p_thr = [p_thr thr(t)];
            %p_d = [p_d d];
            break                                            % Break loop since 10% has been reched
        end
        if round(d,3) < 0.100                                % If density has surpassed 10% (i.e., density is < 10%)
            over = 1;                                        % Indicate that density has been over shot by changing 'over' variable                                 
            t = t-1;                                         % Go back one threshold and re-run (in order to save connectome that has not over-shot 10%)
        end 
    end
end

% Print results
disp(sprintf('At this threshold, a mean density of %g%% is produced across the sample.',100*mean(density)));

