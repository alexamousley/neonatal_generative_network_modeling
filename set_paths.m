%% Set paths to load published data into scripts
%{

Written by Alexa Mousley, MRC Cognition and Brain Sciences Unit
Email: alexa.mousley@mrc-cbu.cam.ac.uk

Edit this script to contain the paths to the data you would like to use to
run the scripts.

%}

%% Function/Toolbox paths

% Add Brain Connectivity Toolbox
%{
All graphy theory measures were calculated using the Brain Connectivity Toolbox
(https://sites.google.com/site/bctnet/home?authuser=0)

Toolbox publication:
Rubinov, M., & Sporns, O. (2010). Complex network measures of brain 
connectivity: uses and interpretations. Neuroimage, 52(3), 1059-1069.
%}
addpath('/set/path/to/BCT/');  

% Add path to consensus network function
%{
Consensus network function can be found here: 
https://www.brainnetworkslab.com/coderesources

The function is documented in this publication:
Betzel, R. F., Griffa, A., Hagmann, P., & Mišić, B. (2019). 
Distance-dependent consensus thresholds for generating group-representative
structural brain networks. Network neuroscience, 3(2), 475-496.

%}
addpath('/set/path/to/consensus/function/');     

%% Data paths

% Add path to demographic data
addpath('/path/.../demographics'); 
% Add path to all observed networks 
addpath('/path/.../observed_networks'); 
% Add path to all simulated networks 
addpath('/path/.../simulated_networks'); 
% Add path to atlas data 
addpath('/path/.../atlas');     
% Add path to derived data (e.g., graph theory measures)
addpath('/path/.../derived_data'); 
% Add path to density-controlled data
addpath('/path/.../density-controlled_data'); 

% Paths to propensity-matched data
addpath('/path/.../propensity_matched_analysis');              % Networks
addpath('/path/.../propensity_matched_analysis/demographics'); % Demographics
addpath('/path/...propensity_matched_analysis/derived_data'); % Derived data 
