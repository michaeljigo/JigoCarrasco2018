function anal_dPrime(subj,expName,nSets,doFit)
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
data = ['../data/',subj,'/',expName,'/cue/'];

% separate the stim files into sets
if ieNotDefined('nSets') % a set = [exo neutral endo]
    nSets = 1; % group all the data into one set
end
% parse files into cueing conditions
parsedFiles = parseFiles(data,nSets,{'response','parameter',...
    'brokenTrial.trialIdx'});

% update nSets based on how many sets could actually be made
nSets = length(parsedFiles);

% initialize figure
eccFig = figure('Name',[subj,' behavior']);

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
    options.bootstrapIterations = 1000;
    
    % compute hit and false alarm rate
    temp = condParser(resp,target,cue,ecc,options);
    % compute d'
    dPrime.perf = squeeze(diff(norminv(temp.perf),[],1));
    dPrime.err = cellfun(@(h,f) sort(norminv(h)-norminv(f)),...
        squeeze(temp.bootDistribution(2,:,:)),squeeze(temp.bootDistribution(1,:,:)),...
        'UniformOutput',false);
    dPrime.err = cellfun(@(d) d([0.025*options.bootstrapIterations, ...
        (1-0.025)*options.bootstrapIterations]),dPrime.err,'UniformOutput',false);
    dPrime.factorLabels = temp.factorLabels;
    
    %% Plotting
    % d'
    figure(eccFig);
    lineCol = 'kr';
    fitX = linspace(0,7.2,1000);
    ecc = unique(setData.targetEcc);
    for c = 1:length(lineCol)
        errorBounds = cell2mat(dPrime.err(c,:));
        scatter(ecc,dPrime.perf(c,:),50,lineCol(c),'filled'); hold on
        errorbar(ecc,dPrime.perf(c,:),abs(errorBounds(1,:)-dPrime.perf(c,:)),...
            abs(errorBounds(2,:)-dPrime.perf(c,:)),'MarkerSize',1,...
            'Linestyle','none','color',lineCol(c));
        
        % fit data
        if doFit.y
            switch doFit.fittype
                case 'polynomial'
                    polyOrder = 3;
                    grpFit = doFit.grpCoeff;
                    
                    % check if fit should be constrained or not
                    if isempty(grpFit)
                        fitCoeff = polyfit(ecc,dPrime.perf(c,:),polyOrder);
                    else
                        contrasts = eye(polyOrder+1);
                        % adjust contrasts to get inequalities correct (-1 corresponds
                        % to a "greater than" inequality)
                        contrasts = contrasts.*(repmat(grpFit,1,size(contrasts,2)));
                        
                        % set constraints to all be 0 such that the coefficients are
                        % either held below or above 0
                        constraints = zeros(size(contrasts,1),1);
                        
                        % make matrix for lsqlin
                        C = ones(length(ecc),polyOrder+1);
                        for o = polyOrder:-1:1
                            C(:,o) = ecc'.*C(:,o+1);
                        end
                        
                        % fit
                        fitCoeff = lsqlin(C,dPrime.perf(c,:),contrasts,constraints);
                    end
                    fitY = polyval(fitCoeff,fitX);
                    
                    % calculate r2 of model
                    r2(c) = calcR2(dPrime.perf(c,:),polyval(fitCoeff,ecc));
                    
                    % get peak using model
                    [~, peakEcc(c)] = max(fitY);
                    peakEcc(c) = fitX(peakEcc(c));
                case 'spline'                    
                    fitCoeff = fit(ecc',dPrime.perf(c,:)','smoothingspline',...
                        'smoothingParam',doFit.smoothParam);
                    fitY = feval(fitCoeff,fitX);
                    r2(c) = calcR2(dPrime.perf(c,:),feval(fitCoeff,ecc)');
                    % get peak of spline
                    [tempMax, peakEcc(c)] = max(fitY);
                    peakEcc(c) = fitX(peakEcc(c));
            end
            % plot model
            plotLine(c) = plot(fitX,fitY,[lineCol(c),'-']);
            % plot peak
            set(gca,'YLim',[-1 3.5]);
            vline(peakEcc(c),[lineCol(c),'-']);
        else
            plotLine(c) = plot(ecc,dPrime.perf(c,:),[lineCol(c),'-']);
        end
    end
    set(gca,'XLim',[-0.8 8]);
    set(gca,'YLim',[-1 3.5]);
    xlabel('eccentricity (dva)');
    ylabel('d prime');
    legend(plotLine,dPrime.factorLabels.factor2);
end

%% Save result
if nSets==1
    groupFolder = ['../data/group/',expName,'/dPrime/'];
    if ~exist(groupFolder,'dir')
        mkdir(groupFolder);
    end
    
    if doFit.y
        save([groupFolder,subj,'_',doFit.fittype,'.mat'],'dPrime','r2','peakEcc');
    else
        save([groupFolder,subj,'_noFit.mat'],'dPrime');
    end
    fprintf('Saved data.\n');
    
    % save figures
    figFolder = ['../data/figures/dPrime/',expName,'/',subj,'/'];
    if ~exist(figFolder,'dir')
        mkdir(figFolder);
    end
    if doFit.y
        saveas(eccFig,[figFolder,'overall_',doFit.fittype,'.tif']);
    else
        saveas(eccFig,[figFolder,'overall_noFit.tif']);
    end
end