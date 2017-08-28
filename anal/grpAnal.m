% Concatenates and plots the average performance across subjects

function stats = grpAnal(expName,exclude,figExt,doFit)
%% Load subject data
% where subject results are stored
groupDir = ['../data/group/',expName,'/'];
subj = dir ([groupDir,'*noFit.mat']);

if ieNotDefined('exclude')
    exclude = 0;
end
if ieNotDefined('doFit')
    doFit.y = 0;
    doFit.fittype = 'spline';
elseif doFit.y
    subj = dir([groupDir,'*',doFit.fittype,'.mat']);
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
    
    rawPCorr(s,:,:) = pCorrEcc.perf;
    rawRT(s,:,:) = rtEcc.perf;
    rawInt(s,:,:,:) = pCorrInt.perf;
    rawStimVal(s,:,:) = stimVal.perf;
    rawHemi(s,:,:,:) = pCorrHemi.perf;
    rawFixBreaks(s) = fixBreaks;
    rawCueComb(s,:,:,:) = cueComb.perf;
    if doFit.y
        rawPeak(s,:) = peakEcc;
        rawR2(s,:) = r2;
    end
end

if doFit.y
    % store some stats
    stats.fit.rawPeak = rawPeak;
    stats.fit.subjR2 = rawR2;
end

% save condition labels
cues = pCorrEcc.factorLabels.factor1;
ecc = cellfun(@(x) str2num(x),pCorrEcc.factorLabels.factor2);
clearvars -except raw* ecc cues subjNames expName doFit figExt
%% Do group analyses
% prop. corr: Cue x Eccentricity ANOVA
pCorr.cueXecc = buildrmaov2Matrix(rawPCorr);
stats.pCorr.cueXecc = rmaov2(pCorr.cueXecc,'printOut=0');

% prop. corr: Interval x Cue x Eccentricitiy ANOVA
pCorr.intXcueXecc = buildrmaov3Matrix(rawInt);
stats.pCorr.intXcueXecc = rmaov33(pCorr.intXcueXecc,'printOut=0');

% prop. corr: Hemifield x Cue x Eccentricity ANOVA
pCorr.hemiXcueXecc = buildrmaov3Matrix(rawHemi);
stats.pCorr.hemiXcueXecc = rmaov33(pCorr.hemiXcueXecc,'printOut=0');

% rt: Cue x Eccentricity ANOVA
rt.cueXecc = buildrmaov2Matrix(rawRT);
stats.rt.cueXecc = rmaov2(rt.cueXecc,'printOut=0');

% contrast: Cue x Eccentricity ANOVA
contrast.cueXecc = buildrmaov2Matrix(rawStimVal);
stats.contrast.cueXecc = rmaov2(contrast.cueXecc,'printOut=0');

% average across subjects
grpPCorr.avg = squeeze(mean(rawPCorr,1));
grpRT.avg = squeeze(mean(rawRT,1));
grpInt.avg = squeeze(mean(rawInt,1));
grpHemi.avg = squeeze(mean(rawHemi,1));
grpContrast.avg = squeeze(mean(rawStimVal,1));
if doFit.y
    grpPeak.avg = nanmean(rawPeak);
    grpR2.avg = nanmean(rawR2);
end

% within-subject error (have to use for loop because I can't figure out how
% to use reshape and permute to do this in one line)
% loop through each eccentricity and get the error between neutral and exo
for c = 1:size(rawPCorr,3)
    % proportion correct
    temp = squeeze(rawPCorr(:,:,c));
    grpPCorr.err(:,c) = withinSubjErr(temp,0)';
    
    % rt
    temp = squeeze(rawRT(:,:,c));
    grpRT.err(:,c) = withinSubjErr(temp,0)';
    
    % contrast
    temp = squeeze(rawStimVal(:,:,c));
    grpContrast.err(:,c) = withinSubjErr(temp,0)';
    
    % intervals
    for int = 1:size(grpInt.avg,2)
        temp = squeeze(rawInt(:,int,:,:));
        temp = squeeze(temp(:,:,c));
        grpInt.err(int,:,c) = withinSubjErr(temp,0)';
    end
    
    % hemifields
    for h = 1:size(grpHemi.avg,2)
        temp = squeeze(rawHemi(:,h,:,:));
        temp = squeeze(temp(:,:,c));
        grpHemi.err(h,:,c) = withinSubjErr(temp,0)';
    end
end
if doFit.y
    grpPeak.err = withinSubjErr(rawPeak);
    grpR2.err = withinSubjErr(rawR2);
    
    % compute the 68% CI of the peaks
    nSample = 10000;
    errBounds = [0.16*nSample, (1-0.16)*nSample];
    for c = 1:size(rawPeak,2)
        temp = sort(bootstrp(nSample,@(x) mean(x),rawPeak(:,c)));
        grpPeak.ci(:,c) = temp(errBounds);
        grpPeak.distribution(:,c) = temp;
    end
end

%% Plotting
%% Proportion correct
behavFig = figure('Name',['Performance: ',num2str(size(rawPCorr,1)), ' subjects']);
lineCol = 'kr'; fitX = linspace(0,7.2,10000);
areaCol = [0.5 0.5 0.5; 0.5 0 0];

for i = 1:size(grpPCorr.avg,1)
    scatter(ecc,grpPCorr.avg(i,:),50,lineCol(i),'filled'); hold on
    errorbar(ecc,grpPCorr.avg(i,:),grpPCorr.err(i,:),...
        'Linestyle','none','color',lineCol(i));
    
    if doFit.y
        % fit data
        switch doFit.fittype
            case 'polynomial'
                % do cross-validation to find the most appropriate
                % polynomial order
                polyOrder = 1:5;
                cvErr = LOOCV_polynomial(ecc,grpPCorr.avg(i,:),polyOrder,0);
                % choose the polynomial order with the lowest
                % cross-validated error
                [~,idx] = min(cvErr);
                polyOrder = polyOrder(idx);
                
                % use cross-validated fit to fit full dataset
                fitCoeff = polyfit(ecc,grpPCorr.avg(i,:),polyOrder);
                fitY = polyval(fitCoeff,fitX);
                % compute r2 of fit
                r2 = calcR2(grpPCorr.avg(i,:),polyval(fitCoeff,ecc));
                stats.fit.coeff{i} = fitCoeff;
                
                % get peak using model
%                 syms eqnX
                %                 if polyOrder==3
                %                     f = (fitCoeff(1)*eqnX^3)+(fitCoeff(2)*eqnX^2)+(fitCoeff(3)*eqnX)...
                %                         +fitCoeff(4);
                %                 elseif polyOrder==2
                %                     f = (fitCoeff(1)*eqnX^2)+(fitCoeff(2)*eqnX)+fitCoeff(3);
                %                 elseif polyOrder==1
                %                     f = (fitCoeff(1)*eqnX)+fitCoeff(2);
                %                 end
                %                 f1 = diff(f);
                %                 localMaxMin = solve(f1);
                %                 localMaxMin = unique(real(eval(localMaxMin)));
                %                 % check whether the solved values are maxima or minima
                %                 idx = [];
                %                 for m = 1:length(localMaxMin)
                %                     % find if the closest value to the solution in model
                %                     [~,idx(m)] = min(abs(fitX-localMaxMin(m)));
                %                     tempMax(m) = fitY(idx(m));
                %                 end
                %                 peakEcc(i) = fitX(idx(tempMax==max(tempMax)));
                [tempMax,peakEcc(i)] = max(fitY);
                peakEcc(i) = fitX(peakEcc(i));
                
                % plot peak 95% CI of individual peaks
                forCI = linspace(grpPeak.ci(1,i),grpPeak.ci(2,i),100);
                aPlot(i) = area(forCI,polyval(fitCoeff,forCI));
                set(aPlot(i),'FaceColor',areaCol(i,:));
                set(aPlot(i),'FaceAlpha',0.1);
                set(aPlot(i),'EdgeAlpha',0);
            case 'spline'
                fitCoeff = fit(ecc',grpPCorr.avg(i,:)',...
                    'smoothingspline','smoothingParam',doFit.smoothParam);
                fitY = feval(fitCoeff,fitX);
                r2 = calcR2(grpPCorr.avg(i,:),feval(fitCoeff,ecc)');
%                 get peak of spline
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
        
        % plot peak
%         line([peakEcc(i) peakEcc(i)],[0 max(tempMax)],'color',lineCol(i),...
%             'Linestyle','-');
        
        % calculate r2 of model and store coefficients of fit
        stats.fit.r2(i) = r2;
    else
        plotLine(i) = plot(ecc,grpPCorr.avg(i,:),[lineCol(i),'-']);
    end
end
% set title and pretty up figure
if doFit.y
    title(sprintf('cue=%.3f; ecc=%.3f; cueXecc=%.3f; neutR2=%.2f; cueR2=%.2f',stats.pCorr.cueXecc.iv1{3},...
        stats.pCorr.cueXecc.iv2{3},stats.pCorr.cueXecc.iv12{3},stats.fit.r2(1),stats.fit.r2(2)));
else
    title(sprintf('cue=%.3f; ecc=%.3f; cueXecc=%.3f;',stats.pCorr.cueXecc.iv1{3},...
        stats.pCorr.cueXecc.iv2{3},stats.pCorr.cueXecc.iv12{3}));
end
set(gca,'XLim',[-0.1 8]);
set(gca,'YLim',[0.6 0.9]);
ylabel('proportion correct');
xlabel('eccentricity');
legend(plotLine,cues);

%% RT
rtFig = figure('Name',['RT: ',num2str(size(rawPCorr,1)), ' subjects']);
for i = 1:size(grpRT.avg,1)
    scatter(ecc,grpRT.avg(i,:),20,lineCol(i),'filled'); hold on
    errorbar(ecc,grpRT.avg(i,:),grpRT.err(i,:),'Linestyle','none',...
        'color',lineCol(i));
    plotLine(i) = plot(ecc,grpRT.avg(i,:),[lineCol(i),'-']);
end
% set title and pretty up figure
title(sprintf('cue=%.3f; ecc=%.3f; cueXecc=%.3f',stats.rt.cueXecc.iv1{3},...
    stats.rt.cueXecc.iv2{3},stats.rt.cueXecc.iv12{3}));
set(gca,'XLim',[-0.1 8]);
set(gca,'YLim',[0.15 0.35]);
ylabel('RT (s)');
xlabel('eccentricity');
legend(plotLine,cues);

%% Interval results
intFig = figure('Name','Intervals');
for int = 1:size(grpInt.avg,1)
    subplot(1,2,int);
    for i = 1:size(grpInt.avg,2)
        perf = squeeze(grpInt.avg(int,i,:));
        err = squeeze(grpInt.err(int,i,:));
        scatter(ecc,perf,20,lineCol(i),'filled'); hold on
        errorbar(ecc,perf,err,'Linestyle','none','color',lineCol(i));
        plot(ecc,perf,[lineCol(i),'-']);
    end
    title(sprintf('Interval #%i',int));
    set(gca,'XLim',[-0.1 8]);
    set(gca,'YLim',[0.5 1]);
    ylabel('proportion correct');
    xlabel('eccentricity');
end
legend(plotLine,cues);

% make plots directly comparing neutral and cued performance in each
% interval
lineStyles = {'-' '--'};
directComp = figure('Name','Neutral/cued interval comparison');
for cue = 1:size(grpInt.avg,2)
    subplot(1,size(grpInt.avg,2),cue);
    cueData = squeeze(grpInt.avg(:,cue,:));
    for int = 1:size(grpInt.avg,1)
        scatter(ecc,cueData(int,:),20,lineCol(cue),'filled'); hold on
        plotLine(int) = plot(ecc,cueData(int,:),[lineCol(cue),lineStyles{int}]);
        errorbar(ecc,cueData(int,:),squeeze(grpInt.err(int,cue,:)),...
            'Linestyle','none','color',lineCol(cue));
    end
    set(gca,'YLim',[0.5 1]);
    set(gca,'XLim',[-0.1 8]);
    ylabel('proportion correct');
    xlabel('eccentricity');
    title(cues{cue});
    legend(plotLine,{'Interval #1' 'Interval #2'});
end
% %% Contrast results
% contrastFig = figure('Name','Contrast');
% for i = 1:size(grpContrast.avg,1)
%     scatter(ecc,grpContrast.avg(i,:),20,lineCol(i),'filled'); hold on
%     errorbar(ecc,grpContrast.avg(i,:),grpRT.err(i,:),'Linestyle','none',...
%         'color',lineCol(i));
%     plotLine(i) = plot(ecc,grpContrast.avg(i,:),[lineCol(i),'-']);
% end
% % set title and pretty up figure
% title(sprintf('cue=%.3f; ecc=%.3f; cueXecc=%.3f',stats.contrast.cueXecc.iv1{3},...
%     stats.contrast.cueXecc.iv2{3},stats.contrast.cueXecc.iv12{3}));
% set(gca,'XLim',[-0.1 8]);
% set(gca,'YLim',[0 1]);
% ylabel('gabor contrast');
% xlabel('eccentricity');
% legend(plotLine,cues);

% %% Hemifield results
% hemiFig = figure('Name','Hemifields');
% hemi = {'left hemi' 'right hemi'};
% for h = 1:size(grpHemi.avg,1)
%     subplot(1,2,h);
%     for i = 1:size(grpHemi.avg,2)
%         perf = squeeze(grpHemi.avg(h,i,:));
%         err = squeeze(grpHemi.err(h,i,:));
%         scatter(ecc,perf,20,lineCol(i),'filled'); hold on
%         errorbar(ecc,perf,err,'Linestyle','none','color',lineCol(i));
%         %         [~,~,out] = doRegression(ecc,perf,2);
%         %         plot(ecc,out,[lineCol(i),'-']);
%         plot(ecc,perf,[lineCol(i),'-']);
%     end
%     title(sprintf('%s',hemi{h}));
%     set(gca,'XLim',[-0.1 8]);
%     set(gca,'YLim',[0.5 1]);
%     ylabel('proportion correct');
%     xlabel('eccentricity');
% end
% legend(plotLine,cues);

%% Scatter plot of neutral vs. cued performance for individual subjects
% Set colors eccentricities (darker colors=fovea; lighter=periphery)
% eccColors = linspace(0,0.8,7);
eccColors = [255 0 0; 255 0.4*255 0.3*255;8 69 148;33 113 181;66 146 198;...
    107 174 214;158 202 225]/255;
subjCueEff = figure('Name','Subject-level cueing effects');

for s = 1:size(rawPCorr,1)
    for e = 1:size(rawPCorr,3)
        forLeg(e) = scatter(rawPCorr(s,1,e),rawPCorr(s,2,e),40,eccColors(e,:),...
            'filled'); hold on
    end
end
h = plot([-1 1],[-1 1],'k--'); hold on
set(h,'Linewidth', 1.1);
set(gca,'YLim',[0.4 1]);
set(gca,'XLim',[0.4 1]);
axis square
legend(forLeg,cellfun(@(x) num2str(x),num2cell(ecc),'UniformOutput',false));

xlabel('NEUTRAL performance');
ylabel([upper(cues{2}), ' performance']);

%% Plot peak calculation across subjects
if doFit.y
    % peaks
    [~,stats.peak.p,~,stats.peak.info] = ttest(rawPeak(:,1),rawPeak(:,2));
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
figFolder = ['../data/figures/group/',expName,'/'];
if ~exist(figFolder,'dir')
    mkdir(figFolder);
end
if ieNotDefined('figExt')
    figExt = [];
end
saveas(intFig,[figFolder,'interval',figExt,'.tif']);
saveas(rtFig,[figFolder,'rt',figExt,'.tif']);
saveas(subjCueEff,[figFolder,'subjCueEff',figExt,'.tif']);
saveas(directComp,[figFolder,'intervalComparison',figExt,'.tif']);
if doFit.y
    saveas(behavFig,[figFolder,'overall_',doFit.fittype,figExt,'.tif']);
    saveas(r2Fig,[figFolder,'subjR2_',doFit.fittype,figExt,'.tif']);
    saveas(peakFig,[figFolder,'peaks_',doFit.fittype,figExt,'.tif']);
else
    saveas(behavFig,[figFolder,'overall',figExt,'.tif']);
end
