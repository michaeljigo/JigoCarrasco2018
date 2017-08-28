% Concatenates and plots the average performance across subjects

function stats = grpAnal_dPrime(expName,exclude,figExt,doFit)
%% Load subject data
% where subject results are stored
groupDir = ['../data/group/',expName,'/dPrime/'];
subj = dir([groupDir,'*noFit.mat']);

if ieNotDefined('exclude')
    exclude = 0;
end
if ieNotDefined('doFit')
    doFit.y = 0;
    doFit.fittype = 'spline';
elseif doFit.y
    subj = dir([groupDir,'*',doFit.fittype,'.mat']);
    peakDir = [groupDir,'bootStrap/'];
end

% get subject initials to make sure correct subjects are being used in
% group analysis
subjNames = arrayfun(@(x) x.name(1:strfind(x.name,'.')-1), subj,...
    'UniformOutput',false);
% display subject names
disp('%%%%% Subjects %%%%%');
cellfun(@(x) fprintf('%s \n',x),subjNames);
if exclude
    try
        subjExc = eval(input(sprintf('Enter the subjects you want excluded (use cell array of strings): '),'s'));
        subjExc = ismember(subjNames,subjExc);
        subj(subjExc) = [];
        subjNames(subjExc) = [];
    catch
        fprintf('Not removing any subjects.\n');
    end
end

% loop through subjects and load data
for s = 1:length(subj)
    load([groupDir,subj(s).name]);
    rawPerf(s,:,:) = dPrime.perf;
    if doFit.y
        rawPeak(s,:) = peakEcc;
        rawR2(s,:) = r2;
        
        % load the peak distributions for individual subjects
        load([peakDir,subj(s).name]);
        % draw n samples from each subject for averaging across subjects
        n = 10000;
        rawDist(s,1,:) = peakDist.neut(randi(length(peakDist.neut),1,n));
        rawDist(s,2,:) = peakDist.cue(randi(length(peakDist.cue),1,n));
    end
end

if doFit.y
    % store some stats
    stats.fit.rawPeak = rawPeak;
    stats.fit.subjR2 = rawR2;
end

% save condition labels
cues = dPrime.factorLabels.factor2;
ecc = cellfun(@(x) str2num(x),dPrime.factorLabels.factor3);
clearvars -except raw* ecc cues subjNames expName doFit figExt stats

%% Do group analyses
% d': Cue x Eccentricity ANOVA
dPrime.cueXecc = buildrmaov2Matrix(rawPerf);
stats.dPrime.cueXecc = rmaov2(dPrime.cueXecc,'printOut=0');

% average across subjects
grpDPrime.avg = squeeze(mean(rawPerf,1));
if doFit.y
    grpPeak.avg = nanmean(rawPeak);
    grpR2.avg = nanmean(rawR2);
    grpPeak.err = withinSubjErr(rawPeak);
    grpR2.err = withinSubjErr(rawR2);
    
    % compute the 95% CI (i.e., 1 std of a normal distribution) of the peaks
    nSample = 10000;
    errBounds = [0.025*nSample, (1-0.025)*nSample];
    
    % make the group-level distribution for each cueing condition
    temp = sort(squeeze(mean(rawDist,1)),2);
    grpPeak.ci = temp(:,errBounds)'; % each column is a cue condition
end

% within-subject error (have to use for loop because I can't figure out how
% to use reshape and permute to do this in one line)
% loop through each eccentricity and get the error between neutral and exo
for c = 1:size(rawPerf,3)
    % proportion correct
    temp = squeeze(rawPerf(:,:,c));
    grpDPrime.err(:,c) = withinSubjErr(temp,0)';
end

%% Plotting
%% Proportion correct
behavFig = figure('Name',['Performance: ',num2str(size(rawPerf,1)), ' subjects']);
lineCol = 'kr'; fitX = linspace(0,7.2,10000);
areaCol = [0.5 0.5 0.5; 0.5 0 0];
for i = 1:size(grpDPrime.avg,1)
    scatter(ecc,grpDPrime.avg(i,:),50,lineCol(i),'filled'); hold on
    errorbar(ecc,grpDPrime.avg(i,:),grpDPrime.err(i,:),...
        'Linestyle','none','color',lineCol(i));
    
    if doFit.y
        % fit data
        switch doFit.fittype
            case 'polynomial'
                % do cross-validation to find the most appropriate
                % polynomial order
                polyOrder = 3;
                cvErr = LOOCV_polynomial(ecc,grpDPrime.avg(i,:),polyOrder,0);
                % choose the polynomial order with the lowest
                % cross-validated error
                [~,idx] = min(cvErr);
                polyOrder = polyOrder(idx);
                
                % use cross-validated fit to fit full dataset
                fitCoeff = polyfit(ecc,grpDPrime.avg(i,:),polyOrder);
                fitY = polyval(fitCoeff,fitX);
                % compute r2 of fit
                r2 = calcR2(grpDPrime.avg(i,:),polyval(fitCoeff,ecc)); % subtracting 1 for constant term
                stats.fit.coeff{i} = fitCoeff;
                
                % plot within-subject standard error of peaks
                forCI = linspace(grpPeak.avg(i)-grpPeak.err(i),grpPeak.avg(i)+grpPeak.err(i),100);
                aPlot(i) = area(forCI,polyval(fitCoeff,forCI));
                set(aPlot(i),'FaceColor',areaCol(i,:));
                set(aPlot(i),'FaceAlpha',0.1);
                set(aPlot(i),'EdgeAlpha',0);
            case 'spline'
                fitCoeff = fit(ecc',grpDPrime.avg(i,:)',...
                    'smoothingspline','smoothingParam',doFit.smoothParam);
                fitY = feval(fitCoeff,fitX);
                
                r2 = calcR2(grpDPrime.avg(i,:),feval(fitCoeff,ecc)');
                % get peak of spline
                [tempMax, peakEcc(i)] = max(fitY);
                peakEcc(i) = fitX(peakEcc(i));
                
                % plot peak 95% CI of individual peaks
                forCI = linspace(grpPeak.ci(1,i),grpPeak.ci(2,i),100);
                aPlot(i) = area(forCI,feval(fitCoeff,forCI));
                set(aPlot(i),'FaceColor',areaCol(i,:));
                set(aPlot(i),'FaceAlpha',0.1);
                set(aPlot(i),'EdgeAlpha',0);
        end
        % plot fit
        plotLine(i) = plot(fitX,fitY,[lineCol(i),'-']);
        set(plotLine(i),'Linewidth',2);
        
        % calculate r2 of model and store coefficients of fit
        stats.fit.r2(i) = r2;
    else
        plotLine(i) = plot(ecc,grpDPrime.avg(i,:),[lineCol(i),'-']);
    end
end
% set title and pretty up figure
if doFit.y
    title(sprintf('cue=%.3f; ecc=%.3f; cueXecc=%.3f; neutR2=%.2f; cueR2=%.2f',stats.dPrime.cueXecc.iv1{3},...
        stats.dPrime.cueXecc.iv2{3},stats.dPrime.cueXecc.iv12{3},stats.fit.r2(1),stats.fit.r2(2)));
else
    title(sprintf('cue=%.3f; ecc=%.3f; cueXecc=%.3f;',stats.dPrime.cueXecc.iv1{3},...
        stats.dPrime.cueXecc.iv2{3},stats.dPrime.cueXecc.iv12{3}));
end
set(gca,'XLim',[-0.3 8]);
set(gca,'YLim',[0.4 2.2]);
set(gca,'XTick',0:2:8);
set(gca,'YTick',0.5:0.5:2);
set(gca,'TickDir','out');
ylabel('d prime');
xlabel('eccentricity');
legend(plotLine,cues);

% %% Scatter plot of neutral vs. cued performance for individual subjects
% % Set colors eccentricities (darker colors=fovea; lighter=periphery)
% % eccColors = linspace(0,0.8,7);
% eccColors = [255 0 0; 255 0.4*255 0.3*255;8 69 148;33 113 181;66 146 198;...
%     107 174 214;158 202 225]/255;
% subjCueEff = figure('Name','Subject-level cueing effects');
%
% for s = 1:size(rawPerf,1)
%     for e = 1:size(rawPerf,3)
%         forLeg(e) = scatter(rawPerf(s,1,e),rawPerf(s,2,e),40,eccColors(e,:),...
%             'filled'); hold on
%     end
% end
% h = plot([-1 3],[-1 3],'k--'); hold on
% set(h,'Linewidth', 1.1);
% set(gca,'YLim',[-1 3]);
% set(gca,'XLim',[-1 3]);
% axis square
% legend(forLeg,cellfun(@(x) num2str(x),num2cell(ecc),'UniformOutput',false));
%
% xlabel('NEUTRAL performance');
% ylabel([upper(cues{2}), ' performance']);

%% Plot peak calculation across subjects
if doFit.y
    % peaks
    [~,stats.peak.p,~,stats.peak.info] = ttest(rawPeak(:,1),rawPeak(:,2));
    % compute cohen's d
    stats.peak.info.cohenD = abs(stats.peak.info.tstat)/sqrt(size(rawPeak,1));
    peakFig = figure('Name','Subject-level peaks');
    bars(grpPeak.avg,grpPeak.err,[],[],'wr');
    set(gca,'XTickLabel',cues);
    ylabel('Peak eccentricity (dva)');
    set(gca,'YLim',[0 4.5]);
    title(sprintf('p=%.3f',stats.peak.p));
    
    % r2
    r2Fig = figure('Name','Subject-level r2');
    bars(grpR2.avg,grpR2.err,[],[],'wr');
    set(gca,'XTickLabel',cues);
    ylabel('R^2');
    set(gca,'YLim',[0 1]);
end

%% Save figures
figFolder = ['../data/figures/group/dPrime/',expName,'/'];
if ~exist(figFolder,'dir')
    mkdir(figFolder);
end
if ieNotDefined('figExt')
    figExt = [];
end
% saveas(subjCueEff,[figFolder,'subjCueEff',figExt,'.tif']);
if doFit.y
    saveas(behavFig,[figFolder,'overall_',doFit.fittype,figExt,'.pdf']);
    saveas(r2Fig,[figFolder,'subjR2_',doFit.fittype,figExt,'.pdf']);
    saveas(peakFig,[figFolder,'peaks_',doFit.fittype,figExt,'.pdf']);
else
    saveas(behavFig,[figFolder,'overall',figExt,'.pdf']);
end
