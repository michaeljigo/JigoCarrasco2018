% Behavioral data has shown that false alarms increase when the response
% cue points to very eccentric locations. This function will check if this
% increase is due to some weird quirk about the stimuli used or if it is
% due to the subject.

% In essence, I will look at all false alarms, get the noise images that
% were presented, and plot a histogram of the carrier noise that most
% contributed to the false alarms at each eccentricity.

function checkForBadCarrier(subj,expName)

%% Initialization stuff
% set data directory
data = ['../data/',subj,'/',expName,'/cue/'];

% parse files into cueing conditions
setData = parseFiles(data,1,{'response','parameter',...
    'brokenTrial.trialIdx','reactionTime' 'randVars'});

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

% find the conditions present
cueTypes = strfind(expName,'_');
cueTypes = expName(cueTypes(end)+1:end);
possibleConds.names = {'neutral' [cueTypes, ' cue']};
possibleConds.val = [0 1]; % these are a 1-to-1 match with the cond names
cueTypes = unique(setData.cue);

% set the order of condTypes to match the order outputted by condParser
[~,order] = ismember(cueTypes,possibleConds.val);
cueTypes = possibleConds.names(order);

%% condParser stuff
yesResp = logical(setData.response==1);
yesTarg = logical(setData.targetPresent==1);

% false alarms
fa.val = yesResp(~yesTarg);
fa.label = {'no' 'yes'};

% eccentricity condition
ecc = unique(setData.targetEcc);
eccLabel = num2cell(ecc);
eccLabel = cellfun(@(x) num2str(x),eccLabel,'UniformOutput',false);
noTargEcc.val = setData.targetEcc(~yesTarg);
noTargEcc.label = eccLabel;

% cue condition
noTarg = setData.cue(~yesTarg);
noTargCue.val = noTarg;
noTargCue.label = cueTypes;

% carriers
car = setData.carrierFlip(~yesTarg);

carOut = condParser(car,fa,noTargCue,noTargEcc);

%% Make histograms
ncols = length(noTargEcc.label);
nrows = length(noTargCue.label);

for r = 1:nrows
    for c = 1:ncols
        subplot(nrows,ncols,(r-1)*ncols+c);
        pDist(r,c,:) = hist(carOut.raw{2,r,c},1:length(unique(setData.carrierFlip)));
        pDist(r,c,:) = pDist(r,c,:)./sum(pDist(r,c,:));
        hist(carOut.raw{2,r,c},1:length(unique(setData.carrierFlip)));
        set(gca,'YLim',[0 10]); set(gca,'XLim',[0 13]);
        title([noTargCue.label{r},'_',noTargEcc.label{c}],'Interpreter','none');
    end
end
pDist(pDist==0) = 0.1;

% plot the normalized change
figure('Name',subj);
for c = 1:ncols
    normChange = squeeze(pDist(2,c,:)./pDist(1,c,:));
    plot(normChange,[getcolor(c),'-']); hold on
end
legend(noTargEcc.label);
xlabel('carriers');
ylabel('cue/neutral change');

% plot overall distribution for each cue and the proportion changes between
% cue conditions
% figure('Name',subj);
% nrows = 1; ncols = length(noTargCue.label);

% for i = 1:ncols
%     subplot(nrows,ncols,i);
%     pDist(i,:) = hist(cell2mat(squeeze(carOut.raw(2,i,:))'),1:12);
% change pDist to proportion
%     pDist(i,:) = pDist(i,:)/sum(pDist(i,:));
% actually plot histogram
%     hist(cell2mat(squeeze(carOut.raw(2,i,:))'),1:12);
%     set(gca,'XLim',[0 13]);
%     yLim(i,:) = get(gca,'YLim');
%     title(noTargCue.label{i});
% end
% % fix up histograms
% yLim = max(yLim(:));
% for i = 1:ncols
%     subplot(nrows,ncols,i);
%     set(gca,'YLim',[0 yLim]);
% end

% now plot the normalized increase/decrease in carrier presentation between
% cued and neutreal conditions.
% normalized change = cue/neutral.
% figure('Name',subj);
% plot(pDist(2,:)./pDist(1,:),'k*-'); hold on
% xlabel('carriers');
% ylabel('normalized change');