% analyzes the behavior during the endogenous cueing experiment

% INPUT:
% expName = [car(LH)_gab(LH)] -- only one letter in the paranethesis should
% be inputted

% nSets = number of sets, each set will contain each condition type that is
% present in the data directory. For example, if there are three condition
% types, each set will comprise of a trio of conditions.
% Sets will be arranged in stimfile order (ideally this order will
% reflect chronological order but not neccessarily so)

function anal_cue_detection(subj,expName,nSets)

% set data directory
data = ['../data/',subj,'/',expName,'/cue/'];

% separate the stim files into sets
if ieNotDefined('nSets') % a set = [exo neutral endo]
    nSets = 1; % group all the data into one set
end
% parse files into cueing conditions
[parsedFiles x] = parseFiles(data,nSets,{'response','parameter','brokenTrial.trialIdx',...
    'reactionTime'});
% update nSets based on how many sets could actually be made
nSets = length(parsedFiles);

% initialize figure
eccFig = figure('Name',[subj,' behavior']);
hemiFig = figure('Name',[subj, ' hemifields']);
rtFig = figure('Name',[subj, ' rt']);
faFig = figure('Name',[subj, ' false alarms']);
hitRate = figure('Name',[subj ' hit rate']);

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
    yesResp = logical(setData.response==1);
    yesTarg = logical(setData.targetPresent==1);
    
    hit = yesResp(yesTarg);
    fa = yesResp(~yesTarg);
    hitRT = setData.reactionTime(yesTarg);
    
    % eccentricity condition
    ecc = unique(setData.targetEcc);
    eccLabel = num2cell(ecc);
    eccLabel = cellfun(@(x) num2str(x),eccLabel,'UniformOutput',false);
    targEcc.val = setData.targetEcc(yesTarg);
    targEcc.label = eccLabel;
    noTargEcc.val = setData.targetEcc(~yesTarg);
    noTargEcc.label = eccLabel;
    
    % cue condition
    targ = setData.cue(yesTarg);
    noTarg = setData.cue(~yesTarg);
    targCue.val = targ;
    targCue.label = cueTypes;
    noTargCue.val = noTarg;
    noTargCue.label = cueTypes;
    
    % hemifield condition
    targ = setData.targetHemi(yesTarg);
    noTarg = setData.targetHemi(~yesTarg);
    targHemi.val = targ;
    targHemi.label = {'left' 'right'};
    noTargHemi.val = noTarg;
    noTargHemi.label = {'left' 'right'};
    
    % accuracy condition
    acc.val = hit;
    acc.label = {'miss' 'hit'};
    
    %% Eccentricity performance, calculate d' at different eccentricities
    % setup options for bootstrapping
    options.bootstrap = 1;
    options.bootstrapIterations = 100;
    
    % hit rate for different cueing conditions
    hitTargEcc = condParser(hit,targCue,targEcc,options);
    hitTargEcc.perf(hitTargEcc.perf==1) = 0.99;
    
    % calculate false alarm rate for cues
    faCue = condParser(fa,noTargCue,options);
    
    % false alarms as a function of eccentricity
    faCueEcc = condParser(fa,noTargCue,noTargEcc,options);
    
    % Hemifield analysis
    hitTargHemi = condParser(hit,targCue,targHemi,targEcc,options);
    idx1 = hitTargHemi.perf==1; idx0 = hitTargHemi.perf==0;
    idx1Loc = find(idx1); idx0Loc = find(idx0);
    if any([idx1(:); idx0(:)]==1)
        try
            hitTargHemi.perf(idx1) = 1-1/(2*length(hitTargEcc.raw{idx1Loc(1)}));
            hitTargHemi.perf(idx0) = 1/(2*length(hitTargEcc.raw{idx0Loc(1)}));
        catch
        end
    end
    
    % RT analysis
    rtTargEcc = condParser(hitRT,targCue,targEcc,options);
    
    % calculate d' across all eccentricities less than 6 dva for neutral
    % trials to verify that threshold is maintained
    neutralIdx = find(ismember(cueTypes,'neutral'));
    if ~isempty(neutralIdx)
        eccIdx = ecc<3;
        neutDPrime = norminv(mean(cell2mat(hitTargEcc.raw(neutralIdx,eccIdx))))-norminv(faCue.perf(neutralIdx));
        fprintf('Set #%i: Neutral: d prime = %.2f \n',s,neutDPrime);
    end
    
    % calculate confidence intervals using bootstrap samples
    p = 0.025;
    ciBounds = round([p 1-p]*options.bootstrapIterations);
    for c = 1:size(hitTargEcc.raw,1)
        % d'
        tempCI = sort(norminv(cell2mat(hitTargEcc.bootDistribution(c,:)))-...
            norminv(repmat(faCue.bootDistribution{c},1,size(hitTargEcc.raw,2))));
        ci.dPrime{c} = tempCI(ciBounds,:);
        
        % fa
        tempCI = sort(faCue.bootDistribution{c});
        ci.fa(:,c) = tempCI(ciBounds,:);
        tempCI = sort(cell2mat(faCueEcc.bootDistribution(c,:)));
        ci.faEcc{c} = tempCI(ciBounds,:);
        
        
        % rt
        tempCI = sort(cell2mat(rtTargEcc.bootDistribution(c,:)));
        ci.rt{c} = tempCI(ciBounds,:);
        
        % hemifields
        for h = 1:size(hitTargHemi.raw,2)
            tempCI = sort(cell2mat(squeeze(hitTargHemi.bootDistribution(c,c,:))')-...
                repmat(faCue.bootDistribution{c},1,size(hitTargEcc.raw,2)));
            ci.hemi{c,h} = tempCI(ciBounds,:);
        end
    end
    
    %% Plotting
    % plot the performance across eccentricities
    figure(eccFig);
    subplot(nRows,nCols,s);
    lineCol = 'krb';
    % calculate d' for each eccentricity and cue condition
    x = unique(setData.targetEcc);
    dPrime = [];
    for i = 1:size(hitTargEcc.perf,1)
        % plot d'
        dPrime(i,:) = norminv(hitTargEcc.perf(i,:))-norminv(faCue.perf(i));
        plot(x,dPrime(i,:),[lineCol(i),'.']); hold on
        errorbar(x,dPrime(i,:),abs(ci.dPrime{i}(1,:)-dPrime(i,:)),abs(ci.dPrime{i}(2,:)-dPrime(i,:)),'Linestyle','none','color',lineCol(i));
        % fit 2nd polynomial
%         [~,~,fit] = doRegression(x,dPrime(i,:),2);
%         plotLine(i) = plot(x,fit,[lineCol(i),'-']);
        plotLine(i) = plot(x,dPrime(i,:),[lineCol(i),'-']);
    end
    set(gca,'XLim',[-0.1 8]);
    set(gca,'YLim',[-1.5 3]);
    xlabel('eccentricity');
    ylabel('d prime');
    legend(plotLine,cueTypes);
    
    
    % plot RT results
    figure(rtFig);
    subplot(nRows,nCols,s);
    % calculate d' for each eccentricity and cue condition
    x = unique(setData.targetEcc);
    for i = 1:size(hitTargEcc.perf,1)
        % plot avg rt
        y = rtTargEcc.perf(i,:);
        plot(x,y,[lineCol(i),'.']); hold on
        rtLine(i) = plot(x,y,['-',lineCol(i),'.']); hold on
        errorbar(x,y,abs(ci.rt{i}(1,:)-y),abs(ci.rt{i}(2,:)-y),'Linestyle','none',...
            'color',lineCol(i));
        
        % fit 2nd order polynomial
%         [~,~,fit] = doRegression(x,y,2);
%         rtLine(i) = plot(x,fit,['-',lineCol(i)]);

    end
    set(gca,'XLim',[-0.1 8]);
    xlabel('eccentricity');
    ylabel('RT (s)');
    
    
    % plot HEMIFIELD results
    figure(hemiFig);
    % plot each pair curves (left/right hemifield); each cueing condition will
    % be a different color, each hemifield will be a different linestyle
    for i = 1:size(hitTargHemi.perf,2)
        cueData = squeeze(hitTargHemi.perf(:,i,:));
        % Keeping FA constant and observing how the hit rate changes as a
        % function of hemifield. This will point out if it was easier or
        % harder to detect the target in a particular hemifield.
        cueData = norminv(cueData)-repmat(norminv(faCue.perf)',1,size(cueData,2));
        subplot(nSets,length(cueTypes),length(cueTypes)*(s-1)+i);
        for c = 1:size(cueData,1)
            hemiLine(c) = plot(x,cueData(c,:),['.',lineCol(c),'-']); hold on
            errorbar(x,cueData(c,:),abs(ci.hemi{c,i}(1,:)-cueData(c,:)),abs(ci.hemi{c,i}(2,:)-cueData(c,:)),'Linestyle','none',...
            'color',lineCol(c));
        end
        set(gca,'XLim',[-0.1 8]);
        set(gca,'YLim',[-1.5 3]);
        title(hitTargHemi.factorLabels.factor2{i});
        legend(hemiLine,hitTargHemi.factorLabels.factor1);
        xlabel('eccentricity'); ylabel('d prime');
    end
    
    % plot false alarms across eccentricities
    figure(faFig);
    subplot(nRows,nCols,s);
    for i = 1:size(faCueEcc.perf,1)
        y = faCueEcc.perf(i,:);
        plot(x,y,[lineCol(i),'.']); hold on
        errorbar(x,y,abs(ci.faEcc{i}(1,:)-y),abs(ci.faEcc{i}(2,:)-y),'Linestyle','none','color',lineCol(i));
        plotLine(i) = plot(x,y,[lineCol(i),'-']);
    end
    ylabel('false alarm');
    xlabel('eccentricity (dva)');
    set(gca,'YLim',[0 0.5]);
    
    % plot hit rate across eccentricities
    subplot(nRows,nCols,s);
    for i = 1:size(faCueEcc.perf,1)
        figure(hitRate);
        y = hitTargEcc.perf(i,:);
        plot(x,y,[lineCol(i),'.']); hold on
        plotLine(i) = plot(x,y,[lineCol(i),'-']);
    end
    
end
figure(eccFig);
legend(plotLine,targCue.label);
figure(rtFig);
legend(rtLine,targCue.label);

% save results in group folder
if nSets==1
    fa = faCue;
    hit = hitTargEcc;
    rt = rtTargEcc;
    
    groupFolder = ['../data/group/',expName,'/'];
    if ~exist(groupFolder,'dir')
        mkdir(groupFolder);
    end
    save([groupFolder,subj,'.mat'],'fa','hit','rt','faCueEcc','hitTargHemi');
    fprintf('Data saved: %s \n',[groupFolder,subj,'.mat']);
end