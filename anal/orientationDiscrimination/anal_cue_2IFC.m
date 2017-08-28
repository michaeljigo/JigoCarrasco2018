% analyzes the behavior during the endogenous cueing experiment

% INPUT:
% expName = [car(LH)_gab(LH)] -- only one letter in the paranethesis should
% be inputted

% nSets = number of sets, each set will contain each condition type that is
% present in the data directory. For example, if there are three condition
% types, each set will comprise of a trio of conditions.
% Sets will be arranged in stimfile order (ideally this order will
% reflect chronological order but not neccessarily so)

function anal_cue(subj,expName,nSets)

% set data directory
data = ['../data/',subj,'/',expName,'/cue/'];

% separate the stim files into sets
if ieNotDefined('nSets') % a set = [exo neutral endo]
    nSets = 1; % group all the data into one set
end
% parse files into cueing conditions
[parsedFiles x] = parseFiles(data,nSets,{'response','parameter','brokenTrial.trialIdx',...
    'reactionTime' 'randVars'});
% update nSets based on how many sets could actually be made
nSets = length(parsedFiles);

% initialize figure
eccFig = figure('Name',[subj,' behavior']);
hemiFig = figure('Name',[subj, ' hemifields']);
rtFig = figure('Name',[subj, ' rt']);
intFig = figure('Name',[subj, 'intervals']);

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
    
    %% Eccentricity performance, calculate prop. corr at each eccentricities
    % setup options for bootstrapping
    options.bootstrap = 1;
    options.bootstrapIterations = 1000;
    
    % proportion correct for different cueing conditions and eccentricities
    pCorrEcc = condParser(acc,cue,ecc,options);
    pval_pCorr = permutationTest(pCorrEcc.raw,10000,'mean(c2,2)-mean(c1,2)'); % test statistic is exo-neutral
    
    % Hemifield analysis
    pCorrHemi = condParser(acc,cue,hemi,ecc,options);
    
    % RT analysis (only correct trials)
    rtEcc = condParser(setData.reactionTime,cue,ecc,options);
%     allFields = {'perf' 'condLabel' 'raw' 'bootDistribution'};
%     for f = 1:length(allFields)
%         rtEcc.(allFields{f}) = squeeze(rtEcc.(allFields{f})(2,:,:));
%     end
    pval_rt = permutationTest(rtEcc.raw,10000,'mean(c2,2)-mean(c1,2)');
    
    % proportion correct for different cueing conditions split for first
    % and second interval
    pCorrInt = condParser(acc,interval,cue,ecc,options);
    
    % calculate proportion correct across all eccentricities less than
    % 6 dva for neutral trials to verify that threshold (75%) was maintained
    neutralIdx = find(ismember(cueTypes,'neutral'));
    if ~isempty(neutralIdx)
        eccIdx = unique(setData.targetEcc)<3;
        neutPCorr = mean(cell2mat(pCorrEcc.raw(neutralIdx,eccIdx)));
        fprintf('Set #%i: Neutral: pCorr = %.2f \n',s,neutPCorr);
    end
    
    % calculate confidence intervals using bootstrap samples
    p = 0.025;
    ciBounds = round([p 1-p]*options.bootstrapIterations);
    for c = 1:size(pCorrEcc.raw,1)
        % proportion correct
        tempCI = sort(cell2mat(pCorrEcc.bootDistribution(c,:)));
        ci.pCorr{c} = tempCI(ciBounds,:);
        
        % rt
        tempCI = sort(cell2mat(rtEcc.bootDistribution(c,:)));
        ci.rt{c} = tempCI(ciBounds,:);
        
        % hemifields
        for h = 1:size(pCorrHemi.raw,2)
            tempCI = sort(cell2mat(squeeze(pCorrHemi.bootDistribution(c,h,:))'));
            ci.hemi{c,h} = tempCI(ciBounds,:);
        end
        
        % intervals
        for i = 1:size(pCorrInt.raw,1)
            tempCI = sort(cell2mat(squeeze(pCorrInt.bootDistribution(i,c,:))'));
            ci.int{c,i} = tempCI(ciBounds,:);
        end
        
    end
    
    %% Plotting
    % plot the performance across eccentricities
    figure(eccFig);
    subplot(nRows,nCols,s);
    lineCol = 'krb';
    % calculate d' for each eccentricity and cue condition
    x = unique(setData.targetEcc);
    sig = pval_pCorr>0.99|pval_pCorr<0.01;
    for i = 1:size(pCorrEcc.perf,1)
        % plot proportion correct
        scatter(x(sig),pCorrEcc.perf(i,sig),[lineCol(i),'*']); hold on
        scatter(x(~sig),pCorrEcc.perf(i,~sig),[lineCol(i),'.']);
        errorbar(x,pCorrEcc.perf(i,:),abs(ci.pCorr{i}(1,:)-pCorrEcc.perf(i,:)),...
            abs(ci.pCorr{i}(2,:)-pCorrEcc.perf(i,:)),'MarkerSize',1,'Linestyle','none','color',lineCol(i));
        plotLine(i) = plot(x,pCorrEcc.perf(i,:),[lineCol(i),'-']);
    end
    set(gca,'XLim',[-0.1 8]);
    set(gca,'YLim',[0.35 1.05]);
    xlabel('eccentricity');
    ylabel('proportion correct');
    legend(plotLine,cueTypes);
    
    % plot RT results
    figure(rtFig);
    subplot(nRows,nCols,s);
    % calculate d' for each eccentricity and cue condition
    x = unique(setData.targetEcc);
    sig = pval_pCorr>0.99|pval_pCorr<0.01;
    for i = 1:size(rtEcc.perf,1)
        % plot avg rt
        y = rtEcc.perf(i,:);
        scatter(x(sig),y(sig),[lineCol(i),'*']); hold on
        scatter(x(~sig),y(~sig),[lineCol(i),'.']);
        rtLine(i) = plot(x,y,['-',lineCol(i)]); hold on
        errorbar(x,y,abs(ci.rt{i}(1,:)-y),abs(ci.rt{i}(2,:)-y),'Linestyle','none',...
            'color',lineCol(i));
    end
    set(gca,'XLim',[-0.1 8]);
    xlabel('eccentricity');
    ylabel('RT (s)');
    
    % plot HEMIFIELD results
    figure(hemiFig);
    % plot each pair curves (left/right hemifield); each cueing condition will
    % be a different color, each hemifield will be a different linestyle
    for i = 1:size(pCorrHemi.perf,2)
        cueData = squeeze(pCorrHemi.perf(:,i,:));
        subplot(nSets,length(unique(setData.targetHemi)),length(cueTypes)*(s-1)+i);
        for c = 1:size(cueData,1)
            hemiLine(c) = plot(x,cueData(c,:),['.',lineCol(c),'-']); hold on
            errorbar(x,cueData(c,:),abs(ci.hemi{c,i}(1,:)-cueData(c,:)),...
                abs(ci.hemi{c,i}(2,:)-cueData(c,:)),'Linestyle','none',...
                'color',lineCol(c));
        end
        set(gca,'XLim',[-0.1 8]);
        set(gca,'YLim',[0.5 1.05]);
        title(pCorrHemi.factorLabels.factor2{i});
        legend(hemiLine,pCorrHemi.factorLabels.factor1);
        xlabel('eccentricity'); ylabel('proportion correct');
    end
    
    % plot INTERVAL results
    figure(intFig);
    for i = 1:size(pCorrInt.perf,1)
        % plot the results of each interval in separate subplots
        subplot(nSets,2,2*(s-1)+i);
        intData = squeeze(pCorrInt.perf(i,:,:));
        for c = 1:size(intData,1)
            intLine(c) = plot(x,intData(c,:),['.',lineCol(c),'-']); hold on
            errorbar(x,intData(c,:),abs(ci.int{c,i}(1,:)-intData(c,:)),...
                abs(ci.int{c,i}(2,:)-intData(c,:)),'Linestyle','none',...
                'color',lineCol(c));
        end
        set(gca,'XLim',[-0.1 8]);
        set(gca,'YLim',[0.5 1.05]);
        title(pCorrInt.factorLabels.factor1{i});
        legend(hemiLine,pCorrInt.factorLabels.factor2);
        xlabel('eccentricity'); ylabel('proportion correct');
    end
    
end
figure(eccFig);
legend(plotLine,cue.label);
figure(rtFig);
legend(rtLine,cue.label);

% save results in group folder
if nSets==1
    groupFolder = ['../data/group/',expName,'/'];
    if ~exist(groupFolder,'dir')
        mkdir(groupFolder);
    end
%     save([groupFolder,subj,'.mat'],'pCorrEcc','pCorrHemi','rtEcc');
    fprintf('Data saved: %s \n',[groupFolder,subj,'.mat']);
end