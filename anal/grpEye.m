function eye = grpEye(expName)

%1) position during first interval
%2) position during second interval
%3) average position across both intervals
%4) proportion of trials in which either interval had an eye position >1
%5) proportion of trials in which both intervals had an eye position >1

dataDir = ['../data/group/',expName,'/eye/'];
subj = dir([dataDir,'*.mat']);

for s = 1:length(subj)
    load([dataDir,subj(s).name]);
    allInt(s,1,:,:) = int1Pos.perf; %1
    allInt(s,2,:,:) = int2Pos.perf; %2
    allOverall(s,:,:) = overallPos.perf; %3
    allIntDev(s,:) = intervalDeviation; %4
    allTrialDev(s,:) = trialDeviation; %5
end

% do ANOVA
% interval x cue x eccentricity
temp = buildrmaov3Matrix(allInt);
eye.intXcueXecc = rmaov33(temp,'printOut=0');
eye.intervalDeviation = allIntDev;
eye.trialDeviation = allTrialDev;

% plot average trial position
figure('Name','Trial eye position during stimulation')
lineCol = 'kr'; ecc = [0 1.2 2.4 3.6 4.8 6 7.2];
avgOvr = squeeze(mean(allOverall,1));
for c = 1:size(avgOvr,2)
    temp = squeeze(allOverall(:,:,c));
    avgOvrErr(:,c) = withinSubjErr(temp,0)';
end
for c = 1:size(avgOvr,1)
    scatter(ecc,avgOvr(c,:),25,lineCol(c),'filled'); hold on
    errorbar(ecc,avgOvr(c,:),avgOvrErr(c,:),'LineStyle','none','color',lineCol(c));
    plotLine(c) = plot(ecc,avgOvr(c,:),[lineCol(c),'-']);
end
set(gca,'XLim',[-0.1 8]);
set(gca,'YLim',[0 1]);
xlabel('eccentricity');
ylabel('mean distance from fixation during stimulation (dva)');
legend(plotLine,{'neutral' expName(11:end)});
