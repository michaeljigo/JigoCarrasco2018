% usage:   results = condParser(dependentVariable,<varargin>)
% by:      michael jigo
% date:    01/29/17
%          --condParser output now keeps track of the number of elements
%          that were nans (added 03/21/16)
%          --condParser output can now store bootstrap samples (added
%          01/29/17)
%          --condParser now saves the original index of the raw values in
%          each condition (05/04/17)
%          --options.geometricMean=1 will calculate geometric mean and output
%          it in perf.gMean (03/08/18)
% purpose: calculates behavior (e.g., accuracy, reaction time,
%          confidence, etc.) in every combination of inputted conditions.
%
%          If you can see a way to improve this, just let me know :)
% =========================================================================
%          Inputs:
%          - depdentVariable: vector consisting of behavioral responses
%            (e.g. accuracy) on every trial.
%          - varargin: should be at least one CONDITION structure. the OPTIONS
%            structure doesn't have to be included unless you want to specify
%            an option.
%
%          CONDITION
%          The "condition" structure should be set up as follows:
%          CONDITION.val = vector containing the value of the condition/factor
%                          on each trial.
%          CONDITION.label = cell array of strings describing the identity
%                            of each possible value. The order of the
%                            strings should correspond to the order of
%                            unique(CONDITION.val).
%
%          For example, if we have a direction condition with 1=leftward
%          motion and 2=rightward motion, the CONDITION structure will be
%          set up as follows:
%          direction.val = [1 1 1 2 2 1 1];
%          direction.label = {'left' 'right'};
%
%          OPTIONS
%          The "options" structure is an optional input that, if inputted,
%          should have at least one of the following fields:
%          OPTIONS.std: if this field has a value of 1, the standard deviation
%                            for each condition will be calculated.
%          OPTIONS.median: if this field has a value of 1, the median of
%                          behavior will be calculated instead of the mean.
%          OPTIONS.bootstrap: if this field has a value of 1, bootstrap
%                             samples will be made for each condition.
%          OPTIONS.bootstrapIterations: number of bootstrap samples to
%                                       make
%          OPTIONS.bootstrapCI: percentile-based confidence interals
%==========================================================================
%          Output:
%          The results structure will always have a perf, condLabel, and raw
%          field. These fields will contiain:
%          perf: the mean (or median) behavioral performance of the
%                dependentVariable for each condition.
%          condLabel: a string describing the combination of conditions.
%          raw: the raw values in the dependentVariable for each condtion.
%          If options.std==1 then there will be an "std" field
%          containing the standard deviation of the mean for each condition.
%
%          These fields will be organized based on the number of conditions
%          within each condition factor. For example, if there are two
%          factors and factor 1 has 2 conditions and factor 2 has 3
%          conditions then results.perf will be a 2x3 matrix.
%          Alternatively, if there are three factors and factor 1 has 2
%          conditions, factor 2 has 3 conditions, and factor 3 has 4
%          conditions then results.perf will be a 2x3x4 matrix.
%==========================================================================

function ret=condParser(depVar, varargin)

%% Check OPTIONS and CONDITION
% check to see if options was inputted
nInputs = 1:length(varargin); removeCell = 0;
for iInput = 1:length(varargin)
   if ~isfield(varargin{iInput},'val')
      options = varargin{iInput};
      removeCell = iInput;
   end
end
% remove options input from varargin variable
varargin = varargin(nInputs(~ismember(nInputs,removeCell)));

% define possible options
possOpt = {'std' 'median' 'bootstrap' 'bootstrapIterations' 'bootstrapCI' ...
   'geometricMean'};
defaultVal = [0 0 0 10000 0 0];
if ieNotDefined('options')
   % make options structure
   for iCond = 1:length(possOpt)
      options.(possOpt{iCond}) = defaultVal(iCond);
   end
else
   inputOpt = fieldnames(options);
   if ~all(ismember(inputOpt,possOpt))
      %error('check the names of your inputted options')
   end
   % add in defaults for the remaining options
   remainOpt = possOpt(~ismember(possOpt,inputOpt));
   remainVal = defaultVal(~ismember(possOpt,inputOpt));
   for iOpt = 1:length(remainOpt)
      options.(remainOpt{iOpt}) = remainVal(iOpt);
   end
end

if isempty(varargin)
   error('Must provide at least one condition variable');
end

%% Initialize factors (i.e., condition combinations) that will be analyzed
nFactor=length(varargin);
for i=1:nFactor
   cond{i}=varargin{i};
   cond{i}.possVal= unique(cond{i}.val);
   nConditions(i) = length(unique(cond{i}.val));
   possVals{i} = unique(cond{i}.val);
   if length(cond{i}.possVal)~=length(cond{i}.label)
      error(['Condition ', num2str(i), ': label length does not match the number of possible values, check input!']);
   end
   % make fCounter, which stands for factorCounter
   fCounter.(['f',num2str(i)]).idx = 1;
   fCounter.(['f',num2str(i)]).possVal = unique(cond{i}.val);
   fCounter.(['f',num2str(i)]).hasReset = 0;
end

%% Compute factors
nConditions = prod(nConditions); % total number of conditions
condTable = []; factors = fieldnames(fCounter); ret = [];

% make condTable
while ~all(size(condTable)==[nConditions nFactor]) || any(isnan(reshape(condTable,1,numel(condTable))))
   for iFactor = 1:nFactor
      thisRow(iFactor) = fCounter.(factors{iFactor}).possVal(fCounter.(factors{iFactor}).idx);
      condMatrixIdx(iFactor) = fCounter.(factors{iFactor}).idx;
      if iFactor==nFactor
         % calculate output for each combination of conditions
         ret = calculatePerformance(ret,cond,thisRow,condMatrixIdx,depVar,options);
         condTable = [condTable; thisRow];
         thisRow = nan(1,length(factors));
         fCounter.(factors{nFactor}).idx = fCounter.(factors{nFactor}).idx+1;
         if fCounter.(factors{nFactor}).idx>length(fCounter.(factors{nFactor}).possVal)
            fCounter.(factors{nFactor}).idx = 1;
            fCounter.(factors{nFactor}).hasReset = 1;
            for changeFactor = length(factors)-1:-1:1
               if fCounter.(factors{changeFactor+1}).hasReset==1
                  fCounter.(factors{changeFactor}).idx = fCounter.(factors{changeFactor}).idx+1;
                  fCounter.(factors{changeFactor+1}).hasReset = 0;
                  if fCounter.(factors{changeFactor}).idx>length(fCounter.(factors{changeFactor}).possVal)
                     fCounter.(factors{changeFactor}).idx = 1;
                     fCounter.(factors{changeFactor}).hasReset = 1;
                  end
               end
            end
         end
      end
   end
end
return

function perf = calculatePerformance(perf,data,cond,condMatrixIdx,depVar,options)
idxStr = [];
if ieNotDefined('perf')
   perf.perf = [];
   perf.condLabel = {};
   perf.raw = {};
   for iFactor = 1:length(data)
      perf.factorLabels.(['factor',num2str(iFactor)]) = data{iFactor}.label;
   end
   
   % Set options
   if options.std
      perf.std = [];
   end
   if options.bootstrap
      perf.bootstrap = {};
   end
   if options.geometricMean
      perf.gMean = [];
   end
end
condLabel = [];
% determine position of output in matrix
outputPos = '(';
for iIdx = 1:length(condMatrixIdx)
   if iIdx<length(condMatrixIdx)
      outputPos = [outputPos,num2str(condMatrixIdx(iIdx)),','];
   else
      outputPos = [outputPos,num2str(condMatrixIdx(iIdx)),')'];
   end
end
% for cells
cellOutputPos = '{';
for iIdx = 1:length(condMatrixIdx)
   if iIdx<length(condMatrixIdx)
      cellOutputPos = [cellOutputPos,num2str(condMatrixIdx(iIdx)),','];
   else
      cellOutputPos = [cellOutputPos,num2str(condMatrixIdx(iIdx)),'}'];
   end
end
for iFactor = 1:length(cond)
   % make labels
   condLabel = [condLabel,data{iFactor}.label{condMatrixIdx(iFactor)},' / '];
   % make indexing string
   thisStr = ['condIdx{',num2str(iFactor),'}'];
   condIdx{iFactor} = data{iFactor}.val==cond(iFactor);
   if iFactor<length(cond)
      separator = ' & ';
   else
      separator = [];
   end
   idxStr = [idxStr thisStr separator];
end

% remove lagging /
endXIdx = strfind(condLabel,' / ');
condLabel(endXIdx(end):end) = [];

% calculate performance and make variables
evalc('condIdx = eval(idxStr)');
% check if the dependent variable has NaNs
if any(isnan(depVar(condIdx)))
   % if nans are found, store the proportion of trials that were nans
   evalc(['perf.numNaN',eval('outputPos'),' = sum(isnan(depVar(condIdx)))']);
end

if options.median
   evalc(['perf.perf',eval('outputPos'),' = nanmedian(depVar(condIdx))']);
elseif ~options.median
   evalc(['perf.perf',eval('outputPos'),' = nanmean(depVar(condIdx))']);
end

if options.std
   evalc(['perf.std',eval('outputPos'),' = std(depVar(condIdx))']);
end

if options.geometricMean
   evalc(['perf.gMean',eval('outputPos'),' = exp(nanmean(log(depVar(condIdx))))']);
end

if options.bootstrap
   if sum(condIdx)>0
      if numel(depVar(condIdx))==1
         evalc(['perf.bootstrap',eval('cellOutputPos'),' = ',...
            'repmat(depVar(condIdx),[options.bootstrapIterations 1]);']);
      else
         evalc(['perf.bootstrap',eval('cellOutputPos'),' = ',...
            'bootstrp(',num2str(options.bootstrapIterations),',@(x) nanmean(x),depVar(condIdx))']);
      end
   else
      evalc(['perf.bootstrap',eval('cellOutputPos'),' = nan']);
   end
   if options.bootstrapCI>0
      evalc(['temp = sort(perf.bootstrap',eval('cellOutputPos'),')']);
      lowerBound = round((1-options.bootstrapCI)/2*options.bootstrapIterations);
      upperBound = round((options.bootstrapCI+(1-options.bootstrapCI)/...
         2)*options.bootstrapIterations);
      evalc(['perf.bootCI',eval('cellOutputPos'), ' = ', ...
         'temp([lowerBound upperBound],:)']);
   end
end

% Store labels
evalc(['perf.condLabel',eval('cellOutputPos'),' = condLabel']);
evalc(['perf.rawOrgIdx',eval('cellOutputPos'),' = find(condIdx)']);
evalc(['perf.raw',eval('cellOutputPos'),' = depVar(condIdx)']);
