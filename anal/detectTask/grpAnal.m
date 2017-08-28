% Concatenates and plots the average performance across subjects

function grpAnal(expName,exclude)

%% Load subject data
if ieNotDefined('exclude')
    exclude = 0;
end

% where subject results are stored
groupDir = ['../data/group/',expName,'/'];
subj = dir ([groupDir,'*.mat']);

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
    catch
        fprintf('Not removing any subjects.\n');
    end
end

% loop through subjects and load data
for s = 1:length(subj)
    load([groupDir,subj(s).name]);
    rawDPrime(s,:,:) = bsxfun(@minus,norminv(hit.perf),norminv(fa.perf)');
    
    rawFA(s,:) = fa.perf;
    rawRT(s,:,:) = rt.perf;
    rawHemiDPrime(s,1,:,:) = bsxfun(@minus,norminv(squeeze(hitTargHemi.perf(1,:,:))),...
        norminv(fa.perf)');
    rawHemiDPrime(s,2,:,:) = bsxfun(@minus,norminv(squeeze(hitTargHemi.perf(2,:,:))),...
        norminv(fa.perf)');
    rawFAEcc(s,:,:) = faCueEcc.perf;
end
% save condition labels
cues = hit.factorLabels.factor1;
ecc = cellfun(@(x) str2num(x),hit.factorLabels.factor2);
clearvars -except raw* ecc cues

%% Do group analyses
% d': Cue x Eccentricity ANOVA
d.cueXecc = buildrmaov2Matrix(rawDPrime);
d.cueXecc = rmaov2(d.cueXecc,'printOut=0');

% d': Hemifield x Cue x Eccentricity ANOVA
d.hemiXcueXecc = buildrmaov3Matrix(rawHemiDPrime);
d.hemiXcueXecc = rmaov33(d.hemiXcueXecc,'printOut=0');

% rt: Cue x Eccentricity ANOVA
rt.cueXecc = buildrmaov2Matrix(rawRT);
rt.cueXecc = rmaov2(rt.cueXecc,'printOut=0');

% fa: t-test
[~, fa.t, ~, fa.tStats] = ttest(rawFA(:,1), rawFA(:,2));

% fa: Cue x Eccentricity ANOVA
fa.cueXecc = buildrmaov2Matrix(rawFAEcc);
fa.cueXecc = rmaov2(fa.cueXecc,'printOut=0');

% average across subjects
grpDPrime.avg = squeeze(mean(rawDPrime,1));
grpFA.avg = mean(rawFA,1);
grpFAEcc.avg = squeeze(mean(rawFAEcc,1));
grpRT.avg = squeeze(mean(rawRT,1));
grpHemiDPrime.avg = squeeze(mean(rawHemiDPrime,1));

% within-subject error (have to use for loop because I can't figure out how
% to use reshape and permute to do this in one line)
for i = 1:size(rawDPrime,1)
    temp = squeeze(rawDPrime(i,:,:))';
    dPrimeErr(i,:) = temp(:)';
end
% now get error and then transform back into 3D matrix
dPrimeErr = withinSubjErr(dPrimeErr);
dPrimeErr = [dPrimeErr(1:length(dPrimeErr)/2); dPrimeErr(length(dPrimeErr)/2+1:end)];

%% Plotting
% plot behavior
% D PRIME
figure('Name','Sensitivity');
lineCol = 'kr';
for i = 1:size(grpDPrime.avg,1)
    plot(ecc,grpDPrime.avg(i,:),[lineCol(i),'.']); hold on
    %     [~,~,out] = doRegression(ecc,grpDPrime.avg(i,:),2);
    %     plotLine(i) = plot(ecc,out,[lineCol(i),'-']);
    plotLine(i) = plot(ecc,grpDPrime.avg(i,:),[lineCol(i),'-']);
    errorbar(ecc,grpDPrime.avg(i,:),dPrimeErr(i,:),'Linestyle','none','color',lineCol(i));
end
set(gca,'XLim',[-0.1 8]);
set(gca,'YLim',[0 2.5]);
ylabel('d prime');
xlabel('eccentricity');
legend(plotLine,cues);

disp('check');
% RT
% figure('Name','RT');
% for i = 1:size(grpRT,1)
%     plot(ecc,grpRT(i,:),[lineCol(i),'.']); hold on
%     plotLine(i) = plot(ecc,grpRT(i,:),[lineCol(i),'-']);
% end
% set(gca,'XLim',[-0.1 8]);
% ylabel('RT');
% xlabel('eccentricity');
% legend(plotLine,cues);
% 
% % FA
% figure('Name','False alarms');
% h = bar(grpFA);
% set(h,'FaceColor','k');
% set(h,'EdgeColor','k');
% set(gca,'XTickLabel',cues);
% ylabel('false alarm');
