% analyzes the behavior during the endogenous cueing experiment

% INPUT:
% expName = [car(LH)_gab(LH)] -- only one letter in the paranethesis should
% be inputted

% nSets = number of sets, each set will contain each condition type that is
% present in the data directory. For example, if there are three condition
% types, each set will comprise of a trio of conditions.
% Sets will be arranged in stimfile order (ideally this order will
% reflect chronological order but not neccessarily so)

function anal_cue(subj,expName,nSets,isprac,doFit)
if ieNotDefined('isprac')
    % if 1, the function will collapse across all practice trials so that I
    % can compute the number of trials subjects performed prior to the main
    % experimental sessions.
    isprac = 0;
end
if ieNotDefined('doFit')
    % grpFit is a vector populated by 1s or -1s. A 1 in a position
    % constrains the coefficient of that power to less than 0 while a -1
    % forces the coefficient to be greater than 0. For example, a
    % vector of [-1 1 1 1] corresponds to a 3rd-order polynomial where
    % the 3rd-power coefficient (x^3) will be greater than 0 and the other
    % coefficients will be less than 0.
    %     grpFit = [];
    doFit.y = 0;
    doFit.fittype = 'polynomial';
end

% set data directory
if isprac
    data = ['../data/',subj,'/',expName,'/prac/'];
else
    data = ['../data/',subj,'/',expName,'/cue/'];
end

% separate the stim files into sets
if ieNotDefined('nSets') % a set = [exo neutral endo]
    nSets = 1; % group all the data into one set
end
% parse files into cueing conditions
[parsedFiles, stimFiles] = parseFiles(data,nSets,{'response','parameter','brokenTrial.trialIdx',...
    'reactionTime' 'stimVal' 'noise.noiseEcc' 'stair.response'});

% update nSets based on how many sets could actually be made
nSets = length(parsedFiles);

% for practice analysis, just calculate the number of trials for neutral
% and cued trials
if isprac
    prac.neutral = sum(parsedFiles.cue==0);
    prac.cue = sum(parsedFiles.cue==1);
    pracFolder = ['../data/group/',expName,'/prac/'];
    if ~exist(pracFolder,'dir')
        mkdir(pracFolder);
    end
    save([pracFolder,subj,'.mat'],'prac');
    return
end

% because I fucked up in JIS 1st endo session (I made it possible to
% respond during any segment, instead of only after fixation turns green),
% I need to adjust my reaction times to account for the time before the
% final segment
if strcmp(subj,'JIS')
    for s = 1:nSets
        pos = 0;
        for f = 1:length(stimFiles(s).name)
            if ~isempty(strfind(stimFiles(s).name{f},'170614'))
                load([data,stimFiles(s).name{f}]);
                exp = getTaskParameters(myscreen,task);
                for t = 1:length(exp.trials)
                    pos = pos+1;
                    if stimulus.brokenTrial.trialIdx(t)
                        continue
                    else
                        realRT = cumsum(diff(exp.trials(t).segtime));
                        realRT = exp.trials(t).reactionTime(end)-realRT(end);
                        parsedFiles(s).reactionTime(pos) = realRT;
                        parsedFiles(s).response(pos) = exp.trials(t).response(end);
                    end
                end
            end
        end
    end
end

% initialize figure
eccFig = figure('Name',[subj,' behavior']); set(eccFig,'Visible','off');
hemiFig = figure('Name',[subj, ' hemifields']); set(hemiFig,'Visible','off');
rtFig = figure('Name',[subj, ' rt']); set(rtFig,'Visible','off');
intFig = figure('Name',[subj, 'intervals']); set(intFig,'Visible','off');
stairFig = figure('Name',[subj, 'staircase']); set(stairFig,'Visible','off');

nRows = round(sqrt(nSets/1.6));
nCols = ceil(nSets/nRows);

for s = 1:nSets
    setData = parsedFiles(s);
    
    %% Trial removal
    % first, calculate the frequency of fixation breaks
    fixBreaks(s) = mean(setData.trialIdx);
    
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
    
    % noise eccentricity (i.e., the cued non-target eccentricity)
    noiseEcc.val = setData.noiseEcc;
    noiseEcc.label = eccLabel;
    
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
    
    % do median split on contrasts
    [~,idx] = sort(setData.stimVal);
    temp = ones(1,length(setData.stimVal))*2;
    temp(idx(1:round(length(setData.stimVal)/2))) = 1; % lower half
    contrast.val = temp;
    contrast.label = {'lower' 'upper'};
    clear temp
    
    %% Eccentricity performance, calculate prop. corr at each eccentricities
    % setup options for bootstrapping
    options.bootstrap = 1;
    options.bootstrapIterations = 1000;
    options.bootstrapCI = 0.95; % get 95% CI
    
    % get performance for each cue (target and noise) combination
    cueComb = condParser(acc,cue,ecc,noiseEcc);
    
    % proportion correct for different contrast levels
    pCorrContrast = condParser(acc,contrast,cue,ecc);
    
    % stimulus value analysis
    stimVal = condParser(setData.stimVal,cue,ecc);
    
    % proportion correct for different cueing conditions and eccentricities
    pCorrEcc = condParser(acc,cue,ecc,options);
    pCorrEcc.perf(pCorrEcc.perf<0.5) = 0.5;
    pval_pCorr = permutationTest(pCorrEcc.raw,1000,'mean(c2,2)-mean(c1,2)'); % test statistic is exo-neutral
    
    % Hemifield analysis
    pCorrHemi = condParser(acc,hemi,cue,ecc,options);
    
    %     RT analysis (only correct trials)
    rtEcc = condParser(setData.reactionTime,accCond,cue,ecc,options);
    %     only use RTs for correct trials
    theseFields = {'perf' 'condLabel' 'raw' 'bootDistribution' 'rawOrgIdx' ...
        'bootCI'};
    for f = 1:length(theseFields)
        rtEcc.(theseFields{f}) = squeeze(rtEcc.(theseFields{f})(2,:,:));
    end
    pval_rt = permutationTest(rtEcc.raw,1000,'mean(c2,2)-mean(c1,2)');
    
    % proportion correct for different cueing conditions split for first
    % and second interval
    pCorrInt = condParser(acc,interval,cue,ecc,options);
    
    % calculate proportion correct across all eccentricities less than
    % 3 for neutral trials to verify that threshold (75%) was maintained
    neutralIdx = find(ismember(cueTypes,'neutral'));
    if ~isempty(neutralIdx)
        eccIdx = unique(setData.targetEcc)<3;
        neutPCorr = mean(cell2mat(pCorrEcc.raw(neutralIdx,eccIdx)));
        fprintf('Set #%i: Neutral: pCorr = %.2f \n',s,neutPCorr);
    end
    
    %% Plotting
    % plot the performance across eccentricities
    figure(eccFig);
    subplot(nRows,nCols,s);
    lineCol = 'krb';
    % calculate d' for each eccentricity and cue condition
    x = unique(setData.targetEcc); fitX = linspace(0,7.2,100);
    sig = pval_pCorr>0.99|pval_pCorr<0.01;
    for i = 1:size(pCorrEcc.perf,1)
        % get error
        errorBounds = cell2mat(pCorrEcc.bootCI(i,:));
        
        % plot proportion correct
        scatter(x(sig),pCorrEcc.perf(i,sig),50,[lineCol(i),'*']); hold on
        scatter(x(~sig),pCorrEcc.perf(i,~sig),[],lineCol(i),'filled');
        errorbar(x,pCorrEcc.perf(i,:),abs(errorBounds(1,:)-pCorrEcc.perf(i,:)),...
            abs(errorBounds(2,:)-pCorrEcc.perf(i,:)),'MarkerSize',1,'Linestyle','none','color',lineCol(i));
        
        % fit data
        if doFit.y
            switch doFit.fittype
                case 'polynomial'
                    polyOrder = 3;
                    grpFit = doFit.grpCoeff;
                    if isempty(grpFit)
                        fitCoeff = polyfit(x,pCorrEcc.perf(i,:),polyOrder);
                    else
                        contrasts = eye(polyOrder+1);
                        % adjust contrasts to get inequalities correct (-1 corresponds
                        % to a "greater than" inequality)
                        contrasts = contrasts.*(repmat(grpFit,1,size(contrasts,2)));
                        
                        % set constraints to all be 0 such that the coefficients are
                        % either held below or above 0
                        constraints = zeros(size(contrasts,1),1);
                        
                        % make matrix for lsqlin
                        C = ones(length(x),polyOrder+1);
                        for o = polyOrder:-1:1
                            C(:,o) = x'.*C(:,o+1);
                        end
                        
                        % fit
                        fitCoeff = lsqlin(C,pCorrEcc.perf(i,:),contrasts,constraints);
                    end
                    fitY = polyval(fitCoeff,fitX);
                    
                    % calculate r2 of model
                    r2(i) = calcR2(pCorrEcc.perf(i,:),polyval(fitCoeff,x),polyOrder);
                    
                    % get peak using model
%                     syms eqnX
%                     if polyOrder==3
%                         f = (fitCoeff(1)*eqnX^3)+(fitCoeff(2)*eqnX^2)+(fitCoeff(3)*eqnX)...
%                             +fitCoeff(4);
%                     elseif polyOrder==2
%                         f = (fitCoeff(1)*eqnX^2)+(fitCoeff(2)*eqnX)+fitCoeff(3);
%                     end
%                     f1 = diff(f);
%                     localMaxMin = solve(f1);
%                     localMaxMin = unique(round(real(eval(localMaxMin))*10000)/10000);
                    
                    % check whether the solved values are maxima or minima
%                     idx = [];
%                     for m = 1:length(localMaxMin)
%                         % find if the closest value to the solution in model
%                         [~,idx(m)] = min(abs(fitX-localMaxMin(m)));
%                         tempMax(m) = fitY(idx(m));
%                         
%                         % remove duplicate peaks
%                         if range(tempMax)==0
%                             tempMax = unique(tempMax);
%                             idx = unique(idx);
%                         end
%                     end
%                     peakEcc(i) = fitX(idx(tempMax==max(tempMax)));
                    [tempMax, peakEcc(i)] = max(fitY);
                    peakEcc(i) = fitX(peakEcc(i));
                case 'spline'
                    fitCoeff = fit(x',pCorrEcc.perf(i,:)','smoothingspline',...
                        'smoothingParam',doFit.smoothParam);
                    fitY = feval(fitCoeff,fitX);
                    r2(i) = calcR2(pCorrEcc.perf(i,:),feval(fitCoeff,x)');
                    % get peak of spline
                    [tempMax, peakEcc(i)] = max(fitY);
                    peakEcc(i) = fitX(peakEcc(i));
                    
            end
            % plot model
            plotLine(i) = plot(fitX,fitY,[lineCol(i),'-']);
            % plot peak
            set(gca,'YLim',[0.5 1]);
            vline(peakEcc(i),[lineCol(i),'-']);
        else
            plotLine(i) = plot(x,pCorrEcc.perf(i,:),[lineCol(i),'-']);
        end
    end
    hline(0.5,'b--');
    set(gca,'XLim',[-0.1 8]);
    set(gca,'YLim',[0 1]);
    xlabel('eccentricity');
    ylabel('proportion correct');
    legend(plotLine,cueTypes);
    
    % plot RT results
    figure(rtFig);
    subplot(nRows,nCols,s);
    % calculate d' for each eccentricity and cue condition
    x = unique(setData.targetEcc);
    sig = pval_rt>0.99|pval_rt<0.01;
    for i = 1:size(rtEcc.perf,1)
        % get error
        errorBounds = cell2mat(rtEcc.bootCI(i,:));
        
        % plot avg rt
        y = rtEcc.perf(i,:);
        scatter(x(sig),y(sig),50,[lineCol(i),'*']); hold on
        scatter(x(~sig),y(~sig),[],lineCol(i),'filled');
        rtLine(i) = plot(x,y,['-',lineCol(i)]); hold on
        errorbar(x,y,abs(errorBounds(1,:)-y),abs(errorBounds(2,:)-y),...
            'Linestyle','none','color',lineCol(i));
    end
    set(gca,'XLim',[-0.1 8]);
    set(gca,'YLim',[0 1.5]);
    xlabel('eccentricity');
    ylabel('RT (s)');
    
    % plot HEMIFIELD results
    figure(hemiFig);
    for h = 1:size(pCorrHemi.perf,2)
        cueData = squeeze(pCorrHemi.perf(h,:,:));
        cueData(:,1) = pCorrEcc.perf(:,1); % correct zeroEcc performance
        pCorrHemi.perf(h,:,:) = cueData; % update pCorrHemi
        cueData = squeeze(pCorrHemi.perf(h,:,:));
        
        % extract error from condParser output
        cueError = squeeze(cell2mat(pCorrHemi.bootCI(h,:,:)));
        
        subplot(nSets,length(unique(setData.targetHemi)),length(cueTypes)*(s-1)+h);
        for c = 1:size(cueData,1)
            % get error for cue condition
            errorBounds = squeeze(cueError(c,:,:));
            
            scatter(x,cueData(c,:),[],lineCol(c),'filled'); hold on
            hemiLine(c) = plot(x,cueData(c,:),[lineCol(c),'-']);
            errorbar(x,cueData(c,:),abs(errorBounds(1,:)-cueData(c,:)),...
                abs(errorBounds(2,:)-cueData(c,:)),'Linestyle','none',...
                'color',lineCol(c));
        end
        set(gca,'XLim',[-0.1 8]);
        set(gca,'YLim',[0 1]);
        title(pCorrHemi.factorLabels.factor1{h});
        legend(hemiLine,pCorrHemi.factorLabels.factor2);
        xlabel('eccentricity'); ylabel('proportion correct');
    end
    
    % plot INTERVAL results
    figure(intFig);
    for i = 1:size(pCorrInt.perf,1)
        % plot the results of each interval in separate subplots
        subplot(nSets,2,2*(s-1)+i);
        intData = squeeze(pCorrInt.perf(i,:,:));
        
        % get cue errors from condParser output
        cueError = squeeze(cell2mat(pCorrInt.bootCI(i,:,:)));
        for c = 1:size(intData,1)
            % get error for cue condition
            errorBounds = squeeze(cueError(c,:,:));
            
            scatter(x,intData(c,:),[],lineCol(c),'filled'); hold on
            intLine(c) = plot(x,intData(c,:),['.',lineCol(c),'-']);
            errorbar(x,intData(c,:),abs(errorBounds(1,:)-intData(c,:)),...
                abs(errorBounds(2,:)-intData(c,:)),'Linestyle','none',...
                'color',lineCol(c));
        end
        hline(0.5,'b--');
        set(gca,'XLim',[-0.1 8]);
        set(gca,'YLim',[0 1]);
        title(pCorrInt.factorLabels.factor1{i});
        legend(hemiLine,pCorrInt.factorLabels.factor2);
        xlabel('eccentricity'); ylabel('proportion correct');
    end
    
    % plot STAIRCASE
    figure(stairFig);
    subplot(nRows,nCols,s);
    plot(setData.stimVal,'k.-'); hold on
    % make correct responses green dots
    corrIdx = find(acc==1);
    scatter(corrIdx,setData.stimVal(corrIdx),'g','filled');
    set(gca,'YLim',[0.1 1]);
    xlabel('trials'); ylabel('contrast')
    
end
figure(eccFig);
legend(plotLine,cue.label);
figure(rtFig);
legend(rtLine,cue.label);

%% Save results in group folder
if nSets==1
    groupFolder = ['../data/group/',expName,'/'];
    if ~exist(groupFolder,'dir')
        mkdir(groupFolder);
    end
    if doFit.y
        save([groupFolder,subj,'_',doFit.fittype,'.mat'],'pCorrEcc','pCorrHemi','rtEcc','pCorrInt', ...
            'pCorrContrast','stimVal','fixBreaks', 'cueComb','peakEcc','r2');
    else
        save([groupFolder,subj,'_noFit.mat'],'pCorrEcc','pCorrHemi','rtEcc','pCorrInt', ...
            'pCorrContrast','stimVal','fixBreaks', 'cueComb');
    end
    fprintf('Data saved: %s \n',[groupFolder,subj,'.mat']);
    
    % save figures
    figFolder = ['../data/figures/',expName,'/',subj,'/'];
    if ~exist(figFolder,'dir')
        mkdir(figFolder);
    end
    if doFit.y
        saveas(eccFig,[figFolder,'overall_',doFit.fittype,'.tif']);
    else
        saveas(eccFig,[figFolder,'overall_noFit.tif']);
    end
    saveas(hemiFig,[figFolder,'hemifield.tif']);
    saveas(intFig,[figFolder,'interval.tif']);
    saveas(stairFig,[figFolder,'stair.tif']);
    saveas(rtFig,[figFolder,'rt.tif']);
end