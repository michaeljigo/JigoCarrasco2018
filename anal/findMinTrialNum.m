% essentially the same as anal_cue but this function will run bootstrapping
% on subsets of trials to find the minimum number of trials that will give
% similar pattern of results to the complete dataset

function findMinTrialNum(subj,expName)

% set data directory
data = ['../data/',subj,'/',expName,'/cue/'];

% parse files into cueing conditions
setData = parseFiles(data,1,{'response','parameter','brokenTrial.trialIdx',...
    'reactionTime'});

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
% trial accuracy
acc = setData.response==setData.targetInterval;

% eccentricity condition
eccLabel = num2cell(unique(setData.targetEcc));
eccLabel = cellfun(@(x) num2str(x),eccLabel,'UniformOutput',false);
ecc.val = setData.targetEcc;
ecc.label = eccLabel;

% cue condition
cue.val = setData.cue;
cue.label = cueTypes;

% hemifield condition
hemi.val = setData.targetHemi;
hemi.label = {'left' 'right'};

% accuracy condition
accCond.val = acc;
accCond.label = {'incorrect' 'correct'};

% interval condition
interval.val = setData.targetInterval;
interval.label = {'Int1' 'Int2'};

%% Invoke condParser to do its thing

% proportion correct for each cueing condition and eccentricity
options.bootstrap = 1;
options.bootstrapIterations = 1000;
pCorrEcc = condParser(acc,cue,ecc,options);

%% Bootstrap
p = 0.95; % percentile
lb = (1-p)/2; % lower bound of confidence interval
ub = p+(1-p)/2; % upper bound
ciBounds = round([lb ub]*options.bootstrapIterations);

% First bootstrap with all the trials to get the "best" estimate of the CI
for c = 1:size(pCorrEcc.bootDistribution,1)
    temp = sort(cell2mat(pCorrEcc.bootDistribution(c,:)));
    completeCI(c,:,:) = temp(ciBounds,:);
end

% Now bootstrap using subsets of trials that differ in the total number of
% trials
nSubSets = 7;
minTrial = unique(cellfun(@(x) length(x),pCorrEcc.raw))/nSubSets;
maxTrial = minTrial*(nSubSets-1);

% initialize stuff
subCount = 0;
subDist = cell(1,length(minTrial:minTrial:maxTrial));
nTotal = (nSubSets-1)*options.bootstrapIterations; totalCount = 0;
disppercent(-inf);

for n = [1 16 88 96 104] %minTrial:minTrial:maxTrial
    subCount = subCount+1;
    bootPCorr = nan(options.bootstrapIterations,length(cue.label),length(ecc.label));
    for i = 1:options.bootstrapIterations
        % for each number of trials, randomly select the desired subset
        % then calculate the avg. performance for each iteration
        bootPCorr(i,:,:) = cellfun(@(x) mean(x(randi(length(x),[1 n]))),pCorrEcc.raw,...
            'UniformOutput',true);
        totalCount = totalCount+1;
        disppercent(totalCount/nTotal);
    end

    % save distribution
    subDist{subCount} = bootPCorr;
    
    % calculate confidence interval for each number of trials
    bootPCorr = sort(bootPCorr,1);
    subsetCI{subCount}{1} = bootPCorr(ciBounds,1,:);
    subsetCI{subCount}{2} = bootPCorr(ciBounds,2,:);
end
disppercent(inf);


% plot figures
% first plot the complete dataset
colors = 'kr';
figure('Name','complete dataset');
x = unique(ecc.val);
y = pCorrEcc.perf;
for c = 1:size(completeCI,1)
   plotLine(c) = plot(x,y(c,:),[colors(c),'.-']); hold on
   errorbar(x,y(c,:),abs(squeeze(completeCI(c,1,:))'-y(c,:)),...
       abs(squeeze(completeCI(c,2,:))'-y(c,:)),'Linestyle','none',...
       'color',colors(c));
end
xlabel('eccentricity'); ylabel('proportion correct');
set(gca,'yLim',[0.4 1]);
legend(plotLine, {'neutral' 'exo'});

% now plot figures with each number of trials
nrows = 2; ncols = 3;
trials = [1 16 88 96 104]%minTrial:minTrial:maxTrial;
y = cellfun(@(x) squeeze(mean(x,1)),subDist,'UniformOutput',false);
for n = 1:length(subDist)
    if mod(n,nrows*ncols)==1
        figure('Name','trial subsets');
        plotn = 1;
    end
    subplot(nrows,ncols,plotn);
    for c = 1:size(completeCI,1)
        plotLine(c) = plot(x,y{n}(c,:),[colors(c),'.-']); hold on
        lb = squeeze(subsetCI{n}{c});
        ub = lb(2,:);
        lb = lb(1,:);
        errorbar(x,y{n}(c,:),abs(lb-y{n}(c,:)),abs(ub-y{n}(c,:)),...
            'Linestyle','none','color',colors(c));
        set(gca,'XLim',[0 max(x)+0.5]);
        set(gca,'YLim',[0.4 1]);
        xlabel('eccentricity'); ylabel('proportion correct');
    end
    plotn = plotn+1;
    title(sprintf('%i trials',trials(n)));
end
legend(plotLine,{'neutral' 'exo'});
