% analyzes the staircase. computes the contrast at which d'~1.5

% INPUT:
% expName = [car(LH)_gab(LH)] -- only one letter in the paranethesis should
% be inputted

function anal_thresh(subj,expName,threshRun,nSets,varargin)

getArgs(varargin,{'autoSave=0','calcThresh=1'});

if calcThresh
    removeEcc = 6;
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
elseif threshRun==3
    data = ['../data/',subj,'/',expName,'/training/'];
end

% separate the stim files into sets
if ieNotDefined('nSets')
    nSets = 1;
end
[parsedFiles filenames] = parseFiles(data,nSets,{'response','parameter','stimVal',...
    'brokenTrial.trialIdx' 'reactionTime'});
% update nSets with the actual number of sets
nSets = length(parsedFiles);

% initialize figure and analysis results
stair = figure('Name','staircases');
eccFig = figure('Name','d prime @ thresh');
rtFig = figure('Name','rt');
nRows = round(sqrt(nSets/1.6));
nCols = ceil(nSets/nRows);
thresh = nan(1,nSets);
dPrime = thresh;

for s = 1:nSets
    setData = parsedFiles(s);
    
    %% Trial removal
    % remove trials that had their targets presented at eccentricities
    % larger than 6 dva
    removeLargeEcc = setData.targetEcc>=removeEcc & setData.targetPresent==1;
    
    % remove trials where fixation was broken
    if isfield(setData,'trialIdx') % trialIdx is 1 where fixation was broken
        fields = fieldnames(setData);
        fields = setdiff(fields,'trialIdx');
        for f = 1:length(fields)
            setData.(fields{f}) = setData.(fields{f})(setData.trialIdx==0&removeLargeEcc==0);
        end
        % remove the trialIdx field
        setData = rmfield(setData,'trialIdx');
    end
    %% performance analysis
    % make responses into logical index
    setData.response = logical(setData.response==1);
    if calcThresh
        % determine hit rate at the most frequently presented stimulus value;
        % false alarm will be calculated across all no-target trials
        thresh(s) = mode(setData.stimVal);
        threshIdx = setData.stimVal==mode(setData.stimVal);
    else
        threshIdx = logical(ones(1,length(setData.stimVal)));
    end
    
    % index necessary variables used @ threshold (for hit rate)
    threshResp = setData.response(threshIdx);
    threshTargPresent = setData.targetPresent(threshIdx);
    threshTargEcc = setData.targetEcc(threshIdx);
    threshTargPresent = logical(threshTargPresent);
    threshRT = setData.reactionTime(threshIdx);
    
    % set up values for condParser
    hit = threshResp(threshTargPresent);
    fa = setData.response(setData.targetPresent==0);
    targRT = threshRT(threshTargPresent);
    
    % calculate d' across all eccentricities under 6 dva
    dPrime(s) = norminv(mean(hit))-norminv(mean(fa));
    
    % calculate d' at different eccentricities
    targEcc.val = threshTargEcc(threshTargPresent);
    targEcc.label = num2cell(unique(targEcc.val));
    targEcc.label = cellfun(@(x) num2str(x),targEcc.label,'UniformOutput',false);
    
    % plot the performance and RT across eccentricities
    try
        % performance
        hitEcc = condParser(hit,targEcc);
        plotEcc = unique(targEcc.val);
        
        figure(eccFig);
        subplot(nRows,nCols,s);
        % make 1s into 0.99
        hitEcc.perf(hitEcc.perf==1) = 0.99; %%%%%%%%
        y = norminv(hitEcc.perf)-norminv(mean(fa));
        plot(plotEcc,y,'k-'); hold on
        % put number of trials at each eccentricity on plot
        nTrials = cellfun('length',hitEcc.raw);
        for e = 1:length(hitEcc.perf)
            h = text(plotEcc(e),y(e),num2str(nTrials(e)));
            set(h,'HorizontalAlignment','center');
            set(h,'FontSize',14);
        end
        set(gca,'YLim',[-0.5 4]);
        xlabel('eccentricity');
        ylabel('d prime');
        title(['Set #',num2str(s)]);
        
        % rt
        rtEcc = condParser(targRT,targEcc);
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
        xlabel('eccentricity');
        ylabel('reaction time (s)');
        title(['Set #',num2str(s)]);
    catch
        fprintf('Not enough trials to do eccentricity analysis.\n');
    end
    
    %% staircase analysis
    % plot these staircases
    figure(stair);
    subplot(nRows,nCols,s);
    plot(setData.stimVal,'k.-'); hold on
    % make hits green dots
    hitIdx = threshIdx==1 & setData.response==1 & setData.targetPresent==1;
    scatter(find(hitIdx),setData.stimVal(hitIdx),'g','filled');
    new = hline(thresh(s),'g-');
    xlabel('trial #');
    ylabel('gabor contrast');
    title(sprintf('thresh=%.2f  (targPresent=%i; targAbsent=%i)',thresh(s),...
        length(hit),length(fa)));
    set(gca,'YLim',[0 1]);
    fprintf('Set #%i: dPrime=%.1f (hit=%.3f; fa=%.3f)\n',s,dPrime(s),mean(hit),mean(fa));
    
    if threshRun==2
        % plot a horizontal line at the original threshold estimate
        org = hline(finalThresh,'r--');
        title(sprintf('thresh=%.2f; targPresent=%i, targAbsent=%i',thresh(s),length(hit),length(fa)));
        orgD = mean(setData.response(setData.stimVal==finalThresh&setData.targetPresent==1));
        orgD = norminv(orgD)-norminv(mean(fa));
        newD = norminv(mean(hit))-norminv(mean(fa));
        legend([new org],{sprintf('new=%.1f',newD) sprintf('original=%.1f',orgD)});
    else
        legend(new,{'new'});
    end
end

% choose which staircase to save
if ~autoSave
    finalThresh = input('Which staircase to save? ');
    if ~isempty(finalThresh)
        finalThresh = thresh(finalThresh);
        if threshRun==2
            data = orgData;
        end
        save([data,subj,'_thresh.mat'],'finalThresh','dPrime');
        fprintf('saved thresh = %.4f \n',finalThresh);
    else
        fprintf('No threshold saved.\n');
    end
else
    finalThresh = thresh(1);
    save([data,subj,'_thresh.mat'],'finalThresh','dPrime');
    fprintf('saved thresh = %.4f \n',finalThresh);
end