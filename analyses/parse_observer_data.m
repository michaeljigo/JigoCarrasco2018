% Purpose:  Parse an observer's behavioral data and compute behavioral metrics (d-prime, criterion).
%           The resulting central performance drops (CPDs) will be fit with 3rd-order polynomial functions.
%
% By:       Michael Jigo
% Edited:   07.01.21
%
% Input:    subj           subject initials
%           exp_name       'exo' or 'endo'

function parse_observer_data(subj,exp_name)
addpath(genpath('./helperfun'));

   % set data directory
      data = ['../data/',subj,'/',exp_name,'/cue/'];


   % parse raw behavioral files
      parsedFiles = parseFiles(data,1,{'response','parameter','brokenTrial.trialIdx','reactionTime'});


   %% Set up variables for analyses
      % remove trials where fixation was broken
         if isfield(parsedFiles,'trialIdx') % trialIdx is 1 where fixation was broken
            fields = fieldnames(parsedFiles);
            fields = setdiff(fields,'trialIdx');
            for f = 1:length(fields)
               parsedFiles.(fields{f}) = parsedFiles.(fields{f})(parsedFiles.trialIdx==0);
            end
            parsedFiles = rmfield(parsedFiles,'trialIdx');
         end
         parsedFiles.reactionTime(parsedFiles.reactionTime<0) = nan;

      % find the conditions present in parsed structure
         possibleConds.names = {'neutral' [exp_name, ' cue']};
         possibleConds.val = [0 1]; % these are a 1-to-1 match with the cond names
         cueTypes = unique(parsedFiles.cue);

      % set the order of condTypes to match the order outputted by condParser
         [~,order] = ismember(cueTypes,possibleConds.val);
         cueTypes = possibleConds.names(order);

      % set up condParser values 
         % dependent variable: set responses as 1=target present & 0=target absent
            resp = parsedFiles.response;
            resp(parsedFiles.response==2) = 0;

         % target presence (present=1st interval; absetnt=2nd interval)
            target.val = parsedFiles.targetInterval;
            target.val(parsedFiles.targetInterval==2) = 0;
            target.label = {'absent' 'present'};

         % eccentricity
            eccLabel = num2cell(unique(parsedFiles.targetEcc));
            eccLabel = cellfun(@(x) num2str(x),eccLabel,'UniformOutput',false);
            ecc.val = parsedFiles.targetEcc;
            ecc.label = eccLabel;

         % cue condition
            cue.val = parsedFiles.cue;
            cue.label = cueTypes;

   %% Analyze the data
      % compute hit and false alarm rate
         temp = condParser(resp,target,cue,ecc);

      % compute d'
         dPrime.perf = squeeze(diff(norminv(temp.perf),[],1));
         dPrime.perf(dPrime.perf<0) = 0;
         dPrime.factorLabels = temp.factorLabels;

      % compute criterion
         crit = -0.5*(squeeze(norminv(temp.perf(1,:,:))+norminv(temp.perf(2,:,:))));

      % compute reaction time
         %options.geometricMean = 1;
         rt = condParser(log(parsedFiles.reactionTime),cue,ecc);
         rt.perf = exp(rt.perf);


   %% Fit polynomial function to d-prime
      % fitting variables
         model.ecc         = linspace(0,7.2,1e6);
         ecc               = unique(parsedFiles.targetEcc);
         model.polyOrder   = 3; % 3rd-order polynomial will be fit

      % fit to each cueing condition separately
         for c = 1:size(dPrime.perf,1)
            % do fit
               fitCoeff          = polyfit(ecc,dPrime.perf(c,:),model.polyOrder);
               fitY              = polyval(fitCoeff,model.ecc);
               model.dPrime(c,:) = fitY;
   
            % calculate r2 of model
               r2(c)             = calcR2(dPrime.perf(c,:),polyval(fitCoeff,ecc));
   
            % get peak using model
               [~, peakEcc(c)] = max(fitY);
               peakEcc(c) = model.ecc(peakEcc(c));
         end

   %% Save result
      savedir = ['../data/parsed_data/',exp_name,'/'];
      if ~exist(savedir,'dir')
         mkdir(savedir);
      end
      save([savedir,subj,'.mat'],'dPrime','r2','peakEcc','crit','model','ecc','rt');
      fprintf(sprintf('Parsed and saved %s - %s data.\n',subj,exp_name));
