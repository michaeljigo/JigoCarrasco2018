% Purpose:  Wrapper function that run the following analyses for each observer 
%           across both attention conditions (exogenous and endogenous)
%
%           1. Parse observer data                             parse_observer_data
%           2. Bootstrap data for confidence intervals         bootstrap_observer_data
%  
%
% By:       Michael Jigo
% Edited:   07.01.21

function do_analyses
addpath(genpath('./helperfun'));

   % Initialize wrapper variables
      subj     = {'AF' 'AJF' 'CB' 'CP' 'MJ' 'MR' 'KN' 'HL'};
      exp_name = {'exo' 'endo'};

   % Do analyses
   for e = 1:numel(exp_name)
      for s = 1:length(subj)
         % parse observer data
            %parse_observer_data(subj{s},exp_name{e})
      
         % bootstrap data for confidence intervals (this function takes some time to run)
            bootstrap_observer_data(subj{s},exp_name{e})
      end
   end
