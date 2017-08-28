function bootstrapPeaks_dPrime(subj,expName)
% set data directory
data = ['../data/',subj,'/',expName,'/cue/'];

% separate the stim files into sets
if ieNotDefined('nSets') % a set = [exo neutral endo]
    nSets = 1; % group all the data into one set
end
% parse files into cueing conditions
parsedFiles = parseFiles(data,1,{'response','parameter',...
    'brokenTrial.trialIdx'});

% update nSets based on how many sets could actually be made
nSets = length(parsedFiles);

nRows = round(sqrt(nSets/1.6));
nCols = ceil(nSets/nRows);

setData = parsedFiles;

%% Trial removal
% remove trials where fixation was broken
if isfield(setData,'trialIdx') % trialIdx is 1 where fixation was broken
    fields = fieldnames(setData);
    fields = setdiff(fields,'trialIdx');
    for f = 1:length(fields)
        setData.(fields{f}) = setData.(fields{f})(setData.trialIdx==0);
    end
    % remove the trialIdx field
    setData = rmfield(setData,'trialIdx');
end

%% Determine the file types present in setData
% find the conditions present
cueTypes = strfind(expName,'_');
cueTypes = expName(cueTypes(end)+1:end);
possibleConds.names = {'neutral' [cueTypes, ' cue']};
possibleConds.val = [0 1]; % these are a 1-to-1 match with the cond names
cueTypes = unique(setData.cue);

% set the order of condTypes to match the order outputted by condParser
[~,order] = ismember(cueTypes,possibleConds.val);
cueTypes = possibleConds.names(order);

%% Set up values for condParser
% set responses as 1=target present & 0=target absent
resp = setData.response;
resp(setData.response==2) = 0;

% target presence (present=1st interval; absetnt=2nd interval)
target.val = setData.targetInterval;
target.val(setData.targetInterval==2) = 0;
target.label = {'absent' 'present'};

% eccentricity
eccLabel = num2cell(unique(setData.targetEcc));
eccLabel = cellfun(@(x) num2str(x),eccLabel,'UniformOutput',false);
ecc.val = setData.targetEcc;
ecc.label = eccLabel;

% cue condition
cue.val = setData.cue;
cue.label = cueTypes;

%% Analysis
options.bootstrap = 1;
options.bootstrapIterations = 10000;

% compute hit and false alarm rate
temp = condParser(resp,target,cue,ecc,options);
% compute d'
dPrime = squeeze(diff(norminv(temp.perf),[],1));
bootDist = cellfun(@(h,f) sort(norminv(h)-norminv(f)),...
    squeeze(temp.bootDistribution(2,:,:)),squeeze(temp.bootDistribution(1,:,:)),...
    'UniformOutput',false);

% put bottstrapped samples for each eccentricity into cells to easily get 
% the fit for each cue condition
neut = cell2mat(bootDist(1,:));
cue = cell2mat(bootDist(2,:));
neut = mat2cell(neut,ones(1,options.bootstrapIterations),size(dPrime,2));
cue = mat2cell(cue,ones(1,options.bootstrapIterations),size(dPrime,2));

% now fit each row with a third-order polynomial and estimate the peak from
% the fit
polyOrder = 3;
peakDist.neut = cellfun(@(x) getFitPeak(unique(ecc.val),x,polyOrder),neut);
peakDist.cue = cellfun(@(x) getFitPeak(unique(ecc.val),x,polyOrder),cue);

% plot the distributions
figure('Name',subj);
distX = 0:0.05:7.2;
n = hist(peakDist.neut,distX);
c = hist(peakDist.cue,distX);
plot(distX,n,'k-'); hold on
plot(distX,c,'r-');
legend({'neutral' 'cue'});

% save the bootstrapped peak distribution
groupFolder = ['../data/group/',expName,'/dPrime/bootStrap/'];
if ~exist(groupFolder,'dir')
    mkdir(groupFolder);
end

save([groupFolder,subj,'_polynomial.mat'],'peakDist');

fprintf('Saved data.\n');

function peak = getFitPeak(x,y,p)
fitX = linspace(min(x),max(x),10000);
% do fit of x and y with a polynomial with order p
peak = polyfit(x,y,p);
peak = polyval(peak,fitX);
[~,peak] = max(peak);
peak = fitX(peak);
