% Purpose:  Pack up the data for results in Jigo & Carrasco, 2018.
%           This function will create an easy-to-handle structure where data are held.
%
% By:       Michael Jigo
%           07.01.21

function [data params] = pack_experiment_results
addpath(genpath('./helperfun'));

%% Initialize parameters
   params.attntype      = {'exo' 'endo'};          % attention types to process
   params.subjList      = {'AJF' 'AF' 'CB' 'CP' 'MJ' 'MR' 'KN' 'HL'}; 

   
%% Loop through attention types and load respective datasets
   for a = 1:numel(params.attntype)
      %% load each subject's data
      for s = 1:numel(params.subjList)
         dataDir = sprintf('../data/parsed_data/%s/',params.attntype{a});
         load(sprintf('%s%s.mat',dataDir,params.subjList{s}));
      

      %% Experiment parameters
            % eccentricities
            data(a).ecc             = ecc;
   
   
      %% Behavior
         % d-prime
         data(a).dPrime(s,:,:)      = dPrime.perf;

         % cueing effect
         data(a).cueing_effect(s,:) = diff(dPrime.perf);

         % reaction time
         data(a).rt(s,:,:)          = rt.perf;

         % criterion
         data(a).criterion(s,:,:)   = crit;
      

      %% Polynomial fit parameters
         data(a).model.ecc          = model.ecc;
         data(a).model.polyOrder    = model.polyOrder;
         data(a).model.peakEcc(s,:) = peakEcc;

         % load bootstrapped peak eccentricity
            bootdir = sprintf('../data/bootstrap/%s/',params.attntype{a});
            bootval = load([bootdir,params.subjList{s}]);
            data(a).boot.peakEcc(s,1,:)   = bootval.peakDist.neut;
            data(a).boot.peakEcc(s,2,:)   = bootval.peakDist.cue;
      end
   end
