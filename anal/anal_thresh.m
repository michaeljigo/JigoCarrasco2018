% analyzes the staircase. computes the contrast at which proportion correct
% is 75%

% INPUT:
% expName = [car(LH)_gab(LH)] -- only one letter in the paranethesis should
% be inputted

function anal_thresh(subj,expName,threshRun,nSets,varargin)

getArgs(varargin,{'autoSave=0','calcThresh=1'});

if calcThresh
    removeEcc = 3;
else
    removeEcc = inf;
end

% set data directory
if threshRun==1
    data = ['../data/',subj,'/',expName,'/thresh/'];
elseif threshRun==2
    data = ['../data/',subj,'/',expName,'/threshCheck/'];
    % load the final threshold estimate from threshold session
    orgData = ['../data/',subj,'/',expName,'/thresh/'];
    load([orgData,subj,'_thresh.mat']);
end

% separate the stim files into sets
if ieNotDefined('nSets')
    nSets = 1;
end
[parsedFiles, stimFiles] = parseFiles(data,nSets,{'response','parameter','stimVal',...
    'brokenTrial.trialIdx' 'reactionTime'});
% update nSets with the actual number of sets
nSets = length(parsedFiles);

% initialize figure and analysis results
stair = figure('Name','staircases');
eccFig = figure('Name','proportion correct @ thresh');
rtFig = figure('Name','reaction time @ thresh');
intFig = figure('Name','interval analysis');
nRows = round(sqrt(nSets/1.6));
nCols = ceil(nSets/nRows);
thresh = nan(1,nSets);

for s = 1:nSets
    setData = parsedFiles(s);
    
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
    %% performance analysis
    if calcThresh
        % determine hit rate at the most frequently presented stimulus value;
        % false alarm will be calculated across all no-target trials
        thresh(s) = mode(setData.stimVal);
        threshIdx = setData.stimVal==mode(setData.stimVal);
    else
        threshIdx = true(1,length(setData.stimVal));
    end
    
    % index necessary variables used @ threshold (for hit rate)
    threshResp = setData.response(threshIdx);
    threshTargEcc = setData.targetEcc(threshIdx);
    threshTargInterval = setData.targetInterval(threshIdx);
    threshRT = setData.reactionTime(threshIdx);
    
    % calculate proportion correct across trials within the "thresholded"
    % eccentricities
    acc = threshResp==threshTargInterval;
    threshAcc = acc(threshTargEcc<removeEcc);
    pCorr = mean(threshAcc);
    fprintf('Set #%i: proportion correct = %.2f \n',s,pCorr);
    
    % calculate prop corr at different eccentricities
    targEcc.val = threshTargEcc;
    targEcc.label = num2cell(unique(targEcc.val));
    targEcc.label = cellfun(@(x) num2str(x),targEcc.label,'UniformOutput',false);
    
    accCond.val = acc;
    accCond.label = {'incorrect' 'correct'};
    
    targInt.val = threshTargInterval;
    targInt.label = {'1' '2'};
    
    % plot the performance and RT across eccentricities
    try
        % performance
        pCorrEcc = condParser(acc,targEcc);
        plotEcc = unique(targEcc.val);
        
        figure(eccFig);
        subplot(nRows,nCols,s);
        plot(plotEcc,pCorrEcc.perf,'k-'); hold on
        % put number of trials at each eccentricity on plot
        nTrials = cellfun('length',pCorrEcc.raw);
        for e = 1:length(pCorrEcc.perf)
            h = text(plotEcc(e),pCorrEcc.perf(e),num2str(nTrials(e)));
            set(h,'HorizontalAlignment','center');
            set(h,'FontSize',14);
        end
        set(gca,'YLim',[0.5 1.05]);
        xlabel('eccentricity');
        ylabel('proportion correct');
        title(['Set #',num2str(s)]);
        
        % reaction time
        rtEcc = condParser(threshRT,targEcc);
        figure(rtFig);
        subplot(nRows,nCols,s);
        plot(plotEcc,rtEcc.perf,'k-'); hold on
        % put number of trials at each eccentricity on plot
        nTrials = cellfun('length',rtEcc.raw);
        for e = 1:length(rtEcc.perf)
            h = text(plotEcc(e),rtEcc.perf(e),num2str(nTrials(e)));
            set(h,'HorizontalAlignment','center');
            set(h,'FontSize',14);
        end
        set(gca,'YLim',[0 1.5]);
        ylabel('RT (s)');
        xlabel('eccentricity (dva');
        title(['Set #',num2str(s)]);
        
        % interval analysis
        intEcc = condParser(acc,targInt,targEcc);
        figure(intFig);
        subplot(nRows,nCols,s);
        lineCol = 'kb';
        for i = 1:size(intEcc.perf,1)
            plotLine(i) = plot(plotEcc,intEcc.perf(i,:),[lineCol(i),'.-']); hold on
        end
        legend(plotLine,{'1st interval' '2nd interval'});
        set(gca,'YLim',[0.45 1]);
        set(gca,'XLim',[-0.1 8]);
    catch
        fprintf('Not enough trials to do eccentricity analysis.\n');
    end
    
    %% staircase analysis
    % plot these staircases
    figure(stair);
    subplot(nRows,nCols,s);
    plot(setData.stimVal,'k.-'); hold on
    % make "correct" responses green dots
    corrIdx = find(setData.stimVal==thresh(s) & setData.response==setData.targetInterval);
    scatter(corrIdx,setData.stimVal(corrIdx),'g','filled');
    new = hline(thresh(s),'g-');
    xlabel('trial #');
    ylabel('gabor contrast');
    title(sprintf('thresh=%.2f; pCorr=%.2f (# trials=%i)',thresh(s),...
        pCorr,length(acc)));
    set(gca,'YLim',[0 1]);
end

% choose which staircase to save
finalThresh = input('Which staircase to save? ');
if ~isempty(finalThresh)
    % load the posterior of the final stim file to use in the main
    % experiment
    load([data,stimFiles(finalThresh).name{end}]);
    posterior = stimulus.stair.pdf;
    save([data,subj,'_thresh.mat'],'posterior');
    fprintf('saved posterior \n');
else
    fprintf('Nothing saved.\n');
end