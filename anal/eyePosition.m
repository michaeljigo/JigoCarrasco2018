% Check the eye position during each trial (before button response) to
% double-check if my online eyetracking feedback was successful

function eyePosition(subj,expName)

% get files
dataDir = ['../data/',subj,'/',expName,'/cue/'];
files = dir([dataDir,'*.mat']);
stimDist = [];

blink = 0;
sacc = 0;
totalTrials = 0;

%% load both the stimfile (to get trial removal indices) and the eye data
for f = 1:length(files)
    load([dataDir,files(f).name]);
    eye = GozdegetTaskEyeTraces([dataDir,files(f).name]);
    
    %% concatenate eye position and segment times of trials that were unbroken
    for t = 1:eye.nTrials
        % segment times
        segTimes = diff(eye.trials(t).segtime);
        if ~stimulus.brokenTrial.trialIdx(t)
            try
                % get mean position during both fixation periods and
                % stimulation periods
                % fixation #1
                fixTime(:,1) = round(sum(segTimes(1:2))*1000); %start
                fixTime(:,2) = round(sum(segTimes(1:3))*1000); %end
                fixPos = [nanmean(eye.eye.xPos(t,fixTime(1):fixTime(2))), ...
                    nanmean(eye.eye.yPos(t,fixTime(1):fixTime(2)))];
                
                % stimulus #1
                stimTime(:,1) = round(sum(segTimes(1:5))*1000);
                stimTime(:,2) = round(sum(segTimes(1:6))*1000);
                stimPos = [nanmean(eye.eye.xPos(t,stimTime(1):stimTime(2))), ...
                    nanmean(eye.eye.yPos(t,stimTime(1):stimTime(2)))];
                % correct with fixation position
                stimPos = stimPos-fixPos;
                stimDist(end+1,1) = sqrt(sum(stimPos.^2));
                
                % fixatoin #2
                fixTime(:,1) = round(sum(segTimes(1:8))*1000); %start
                fixTime(:,2) = round(sum(segTimes(1:9))*1000); %end
                fixPos = [nanmean(eye.eye.xPos(t,fixTime(1):fixTime(2))), ...
                    nanmean(eye.eye.yPos(t,fixTime(1):fixTime(2)))];
                
                % stimulus #2
                stimTime(:,1) = round(sum(segTimes(1:11))*1000);
                stimTime(:,2) = round(sum(segTimes(1:12))*1000);
                stimPos = [nanmean(eye.eye.xPos(t,stimTime(1):stimTime(2))), ...
                    nanmean(eye.eye.yPos(t,stimTime(1):stimTime(2)))];
                % correct with fixation position
                stimPos = stimPos-fixPos;
                stimDist(end,2) = sqrt(sum(stimPos.^2));
            catch
                stimDist(end+1,:) = nan(1,2);
            end
        else
            % determine if it was a blink break or a saccade break
            
            % segment times
%             startTime = round(sum(segTimes(1:2))*1000);
%             if any(isnan(eye.eye.xPos(t,startTime:end)))
%                 blink = blink+1;
%             else
%                 sacc = sacc+1;
%             end
        end
%         totalTrials = totalTrials+1;
    end
end

%% Separate eye positions as a function of eccentricity, cue, and target interval
% Get condition labels for each trial
setData = parseFiles(dataDir,1,{'parameter','brokenTrial.trialIdx'});

% remove broken trials
if isfield(setData,'trialIdx') % trialIdx is 1 where fixation was broken
    fields = fieldnames(setData);
    fields = setdiff(fields,'trialIdx');
    for f = 1:length(fields)
        setData.(fields{f}) = setData.(fields{f})(setData.trialIdx==0);
    end
    % remove the trialIdx field
    setData = rmfield(setData,'trialIdx');
end

% create condParser variables
cue.val = setData.cue;
cue.label = {'neutral' 'cue'};

eccLabel = num2cell(unique(setData.targetEcc));
eccLabel = cellfun(@(x) num2str(x),eccLabel,'UniformOutput',false);
ecc.val = setData.targetEcc;
ecc.label = eccLabel;

% do analysis separately for interval 1 and interval 2
int1Pos = condParser(stimDist(:,1),cue,ecc);
int2Pos = condParser(stimDist(:,2),cue,ecc);

% analyze the eye position, collapsed across intervals
overallPos = condParser(mean(stimDist,2),cue,ecc);

cuePos = condParser(mean(stimDist,2),cue);

% get the proportion of trials that had eye positions that deviated >1 deg
% from fixation...
% in either interval
thresh = 1;
intervalDeviation = sum(sum(stimDist>thresh,2)>0)/length(stimDist);

% averaged across both intervals in a single trial
trialDeviation = sum(mean(stimDist,2)>thresh)/length(stimDist);

% save the data for group analysis
groupFolder = ['../data/group/',expName,'/eye/'];
if ~exist(groupFolder,'dir')
    mkdir(groupFolder);
end
save([groupFolder,subj,'.mat'],'int1Pos','int2Pos','overallPos',...
    'intervalDeviation', 'trialDeviation','cuePos');
fprintf('Data saved: %s \n',[groupFolder,subj,'.mat']);