% Concatenates and plots the average performance across subjects

function stats = analBothCues_dPrime(exclude,doFit)

%% Load subject data
if ieNotDefined('exclude')
    exclude = 0;
end
if ieNotDefined('doFit')
    doFit.y = 0;
end

% where subject results are stored
expName = {'carL_gabL_exo','carL_gabL_endo'};
cueName = cellfun(@(x) x(11:end),expName,'UniformOutput',false);

for t = 1:length(expName)
    groupDir = ['../data/group/',expName{t},'/dPrime/'];
    subj = dir ([groupDir,'*noFit.mat']);
    
    % get subject initials to make sure correct subjects are being used in
    % group analysis
    subjNames = arrayfun(@(x) x.name(1:strfind(x.name,'.')-1), subj,...
        'UniformOutput',false);
    % display subject names
    disp(['%%%%% ',expName{t},' Subjects %%%%%']);
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
        
        rawDPrime(s,t,:,:) = dPrime.perf;
        rawCueEffect(s,t,:) = diff(dPrime.perf,1);
    end
end

%% Do group analyses
% prop. corr: Task x Cue x Eccentricity ANOVA
temp = buildrmaov3Matrix(rawDPrime);
stats.dPrime.taskXcueXecc = rmaov33(temp,'printOut=0');

%% Calculate error bars for cue effect plot
for c = 1:size(rawCueEffect,3)
    temp = squeeze(rawCueEffect(:,:,c));
    grpCueEffect.err(:,c) = withinSubjErr(temp,0)';
end
grpCueEffect.avg = squeeze(mean(rawCueEffect,1));

%% Plotting
% make scatter plot comparing cueing effect between attention conditions
subjScatter = figure('Name','Subj-level cue-effect');
eccColors = [255 0 0; 255 0.4*255 0.3*255;8 69 148;33 113 181;66 146 198;...
    107 174 214;158 202 225]/255;
cueEff = squeeze(diff(rawDPrime,[],3));
for s = 1:size(cueEff,1)
    for e = 1:size(cueEff,3)
        forLeg(e) = scatter(cueEff(s,1,e),cueEff(s,2,e),50,eccColors(e,:),...
            'filled'); hold on
    end
end
plot([-3 3],[-3 3],'k--');
plot([-2 2],[0 0],'k--');
plot([0 0],[-2 2],'k--');
set(gca,'YLim',[-2 2]); set(gca,'XLim',[-2 2]);
axis square
xlabel('EXO','Interpreter','none');
ylabel('ENDO','Interpreter','none');
title('Cueing effect (valid-neutral)');
legend(forLeg,dPrime.factorLabels.factor3);

% Plot cueing effects at each eccentricity for each cue type
cueEffect = figure('Name','Exo/Endo cueing effect');
ecc = cellfun(@(x) str2num(x),dPrime.factorLabels.factor3);
lineCol = 'kr';
fitX = linspace(0,7.2,100);
for c = 1:size(grpCueEffect.avg,1)
    y = grpCueEffect.avg(c,:);
    yErr = grpCueEffect.err(c,:);
    scatter(ecc,y,50,lineCol(c),'filled'); hold on
    errorbar(ecc,y,yErr,'Linestyle','none','color',lineCol(c));
    
    % FITTING
    if doFit.y
        % fit data
        switch doFit.fittype
            case 'polynomial'
                 % do cross-validation to find the most appropriate
                % polynomial order
                polyOrder = 1;
                cvErr = LOOCV_polynomial(ecc,y,polyOrder,0);
                % choose the polynomial order with the lowest
                % cross-validated error
                [~,idx] = min(cvErr);
                polyOrder = polyOrder(idx);
                
                % use cross-validated fit to fit full dataset
                fitCoeff = polyfit(ecc,y,polyOrder);
                fitY = polyval(fitCoeff,fitX);
                % compute r2 of fit
                r2 = calcR2(y,polyval(fitCoeff,ecc),polyOrder-1);
            case 'spline'
                nTries = 2;
                for n = 1:nTries
                    if n==1
                        fitCoeff = spap2(2,3,ecc,y);
                    else
                        fitCoeff = spap2(newknt(fitCoeff),3,ecc,y);
                    end
                end
                fitY = fnval(fitCoeff,fitX);
                r2 = calcR2(y,fnval(fitCoeff,ecc));
        end
        % plot fit
        plotLine(c) = semilogx(fitX,fitY,[lineCol(c),'-']);
        set(plotLine(c),'Linewidth',2);
        
        % calculate r2 of model and store coefficients of fit
        stats.fit.r2(c) = r2;
        stats.fit.coeff{c} = fitCoeff;
    else
        plotLine(c) = plot(ecc,y,[lineCol(c),'-']);
    end
end

legend(plotLine,cueName);
set(gca,'XLim',[-0.3 8]);
set(gca,'YLim',[-0.4 0.9]);
set(gca,'YTick',-0.3:0.3:0.9);
set(gca,'XTick',0:2:8);
set(gca,'TickDir','out');
xlabel('eccentricity (dva)');
ylabel('cueing effect (valid-neutral)');
% run 2-way ANOVA (cue (Exo vs Endo) X eccentricity)
mat = buildrmaov2Matrix(rawCueEffect);
stats.cueEffect = rmaov2(mat,'printOut=0');
if doFit.y
title(sprintf('exoR2=%.2f; endoR2=%.2f',...
    stats.fit.r2(1),stats.fit.r2(2)));    
else
title(sprintf('cue=%.3f; ecc=%.3f; cueXecc=%.3f',stats.cueEffect.iv1{3},...
    stats.cueEffect.iv2{3},stats.cueEffect.iv12{3}));
end

%% save figures
figFolder = '../data/figures/group/dPrime/';
if ~exist(figFolder,'dir')
    mkdir(figFolder);
end
if ieNotDefined('figExt')
    figExt = [];
end
saveas(subjScatter,[figFolder,'exoVsEndoSubj',figExt,'.pdf']);
if doFit.y
    saveas(cueEffect,[figFolder,'exoVSendoCueEffect_',doFit.fittype,figExt,'.pdf']);
else
    saveas(cueEffect,[figFolder,'exoVSendoCueEffect',figExt,'.pdf']);
end