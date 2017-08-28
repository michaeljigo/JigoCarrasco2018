% analyzes the staircase. computes the contrast at which proportion correct
% is 75%

% INPUT:
% expName = [car(LH)_gab(LH)] -- only one letter in the paranethesis should
% be inputted

function anal_thresh_tilts(subj,expName,threshRun,nSets,varargin)

getArgs(varargin,{'autoSave=0','calcThresh=1'});

% remove eccentricities greater than removeEcc from the performance
% calculation
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
    'brokenTrial.trialIdx' 'reactionTime' 'stair.ccw.pdf'});
% update nSets with the actual number of sets
nSets = length(parsedFiles);

% initialize figure and analysis results
stair = figure('Name','staircases');
eccFig = figure('Name','proportion correct @ thresh');
rtFig = figure('Name','reaction time @ thresh');
nRows = round(sqrt(nSets/1.6));
nCols = ceil(nSets/nRows);

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
    
    % sort stimVal based on staircase type
    vertOrTilt.val = setData.targetTilt;
    vertOrTilt.label = {'vertical' 'tilted'};
    stairType.val = setData.tiltDir;
    stairType.label = {'CCW' 'CW'};
    stairRes = condParser(setData.stimVal,vertOrTilt,stairType);
    % only take tilted trials
    stairRes.raw = stairRes.raw(2,:);
    stairRes.rawOrgIdx = stairRes.rawOrgIdx(2,:);
    % compute mode for each staircase
    thresh(s,:) = cellfun(@(x) mode(x),stairRes.raw,'UniformOutput',false);
    % index trials with the threshold estimate
    if calcThresh
        threshIdx = cellfun(@(x,y,z) z(x==y),stairRes.raw',thresh(s,:)',stairRes.rawOrgIdx','UniformOutput',false);
        stairThresh = inf(1,length(setData.tiltDir));
        % make index for trials at estimated threshold
        for thisThresh = 1:size(threshIdx,1)
            stairThresh(threshIdx{thisThresh}) = thisThresh;
        end
        clear threshIdx
        threshIdx.val = stairThresh;
        threshIdx.label = [stairType.label {'non-threshold'}];
    else
        threshIdx.val = setData.targetTilt;
        threshIdx.label = stairType.label;
    end
    
    % calculate proportion correct across trials within the "thresholded"
    % eccentricities
    acc = setData.response==setData.targetTilt; acc = +acc; % convert to double
    threshAcc = acc;
    threshAcc(setData.targetEcc>removeEcc) = nan;
    threshPerf = condParser(threshAcc,threshIdx);
    fprintf('Set #%i [ccw cw]: [%.2f %.2f] \n',s,threshPerf.perf(1),threshPerf.perf(2));
    
    % calculate prop corr at different eccentricities
    targEcc.val = setData.targetEcc;
    targEcc.label = num2cell(unique(targEcc.val));
    targEcc.label = cellfun(@(x) num2str(x),targEcc.label,'UniformOutput',false);
    
    accCond.val = acc;
    accCond.label = {'incorrect' 'correct'};
    
    % plot the performance and RT across eccentricities
    try
        lineCol = 'kb';
        % performance
        pCorrEcc = condParser(acc,threshIdx,targEcc);
        pCorrEcc.perf = pCorrEcc.perf(1:2,:);
        plotEcc = unique(targEcc.val);
        % get number of trials at each eccentricity
        nTrials = cellfun('length',pCorrEcc.raw);
        
        figure(eccFig);
        subplot(nRows,nCols,s);
        for id = 1:size(pCorrEcc.perf,1)
            perfLine(id) = plot(plotEcc,pCorrEcc.perf(id,:),[lineCol(id),'-']); hold on
            for e = 1:length(pCorrEcc.perf)
                h = text(plotEcc(e),pCorrEcc.perf(id,e),num2str(nTrials(id,e)));
                set(h,'HorizontalAlignment','center');
                set(h,'FontSize',14);
                set(h,'color',lineCol(id));
            end
        end
        set(gca,'YLim',[0.35 1.05]);
        xlabel('eccentricity');
        ylabel('proportion correct');
        title(['Set #',num2str(s)]);
        legend(perfLine,{'CCW' 'CW'});
        
        % reaction time
        rtEcc = condParser(setData.reactionTime,threshIdx,targEcc);
        rtEcc.perf = rtEcc.perf(1:2,:);
        figure(rtFig);
        subplot(nRows,nCols,s);
        for id = 1:size(rtEcc.perf,1)
            scatter(plotEcc,rtEcc.perf(id,:),lineCol(id),'filled'); hold on
            perfLine(id) = plot(plotEcc,rtEcc.perf(id,:),[lineCol(id),'-']);
        end
        set(gca,'YLim',[0 1.5]);
        ylabel('RT (s)');
        xlabel('eccentricity (dva');
        title(['Set #',num2str(s)]);
        legend(perfLine,{'CCW' 'CW'});
    catch
        fprintf('Not enough trials to do eccentricity analysis.\n');
    end
    
    %% staircase analysis
    % plot these staircases
    figure(stair);
    subplot(nRows,nCols,s);
    for id = 1:length(stairRes.perf)
        stairLine(id) = plot(stairRes.raw{id},[lineCol(id),'.-']); hold on
        % make correct responses green dots
        corrIdx = logical(acc(stairRes.rawOrgIdx{id}));
        corrIdx(stairRes.raw{id}~=thresh{s,id}) = 0;
        corrIdx = find(corrIdx);
        scatter(corrIdx,stairRes.raw{id}(corrIdx),'g','filled');
        new(id) = hline(thresh{s,id},[lineCol(id),'--']);
    end
    xlabel('trial #');
    ylabel('gabor tilt (deg)');
    title(sprintf('ccwThresh=%.2f (n=%i); cwThresh=%.2f (n=%i)',...
        thresh{s,1},sum(stairRes.raw{1}==thresh{s,1}),thresh{s,2},...
        sum(stairRes.raw{2}==thresh{s,2})));
    set(gca,'YLim',[0 50]);
    legend(stairLine,{'ccw' 'cw'});
    
    if threshRun==2
        % plot a horizontal line at the original threshold estimate
        for id = 1:length(stairRes.perf)
            org(id) = hline(finalThresh(id),[lineCol(id),'-']);
        end
        legend([new org],{'newCCW' 'newCW' 'oldCCW' 'oldCW'});
    end
end

% choose which staircase to save
if ~autoSave
    finalThresh = input('Which staircase to save? ');
    if ~isempty(finalThresh)
        %         finalThresh = cell2mat(thresh(finalThresh,:));
        if threshRun==2
            data = orgData;
        end
        %         save([data,subj,'_thresh.mat'],'finalThresh');
        %         fprintf('saved thresh [ccw cw]: [%.4f %.4f] \n',finalThresh);
        % load the posterior from the staircase that is going to be saved
        load([data,stimFiles(finalThresh).name{end}]);
        ccwPosterior = stimulus.stair.ccw.pdf;
        cwPosterior = stimulus.stair.cw.pdf;
        save([data,subj,'_thresh.mat'],'ccwPosterior','cwPosterior');
        fprintf('Posterior saved. \n');
    else
        fprintf('No threshold saved.\n');
    end
else
    finalThresh = thresh(1);
    save([data,subj,'_thresh.mat'],'finalThresh');
    fprintf('saved thresh = %.4f \n',finalThresh);
end