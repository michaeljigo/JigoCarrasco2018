% runs exogenous cueing texture segmentation task with second-order
% textures

function runExoTexDetect_withRespCue(myscreen)
global stimulus % initStimulus in main driver function

%% Set parameters for stimuli
% exogenous cue
stimulus.exoCue.rMin = 0.35; % starting point of line
stimulus.exoCue.rMax = stimulus.exoCue.rMin*2; % extent of line
stimulus.exoCue.yDistFromTarg = 1.6;
stimulus.exoCue.size = 3;
stimulus.exoCue.color = [0 0.5 0.13];%*1.05;

% fixation
stimulus.fixation.size = [0.5 0.125];
stimulus.fixation.edgeOffset = 0.07;
stimulus.fixation.color = [0 0 0];
stimulus.fixation.loc = [-0.02 0.02]; % for some reason the fixation cross isn't perfectly centered, this fixes it

% location cue
stimulus.locCue.size = 0.75;
stimulus.locCue.length = 0.7;
stimulus.locCue.color = [0.75 0.75 0.75];

% neutral cue
% make cue just horizontal lines above and below horizontal meridian at
% the same distance as the exo cue. this way the location of the
% stimulus energy is equated
stimulus.frame.x0 = stimulus.carrier.width/2*-1;
stimulus.frame.x1 = stimulus.carrier.width/2;
stimulus.frame.y = stimulus.carrier.height/2;

% this variable lets me run subjects through a quick refresher session
% between eccentricity sets
if ~isfield(stimulus,'targetInterval')
    stimulus.targetInterval = [1 2];
    stimulus.numTrials = 112;
end

%% Initialize stimuVal
stimulus.stimVal = [];

%% Set up task parameters
eyeContingent = stimulus.eye.eyeContingent; % screen to make sure subject is fixating
fixation = 0.5;
cue = 0.04;
cueISI = 0.05;
stim = 0.1; % 10 frames (larger than maximum duration from Yeshurun & Carrasco, 2000)
resp = inf;
respISI = 0.2;
makeTex = inf;
locCue = 0.2;

task{1}.seglen = [makeTex eyeContingent fixation cue cueISI stim respISI locCue ...
    fixation cue cueISI stim respISI locCue respISI resp];
task{1}.getResponse = [zeros(1,length(task{1}.seglen)-1) 1];
task{1}.waitForBacktick = 0;
task{1}.numTrials = 56;%stimulus.numTrials; % 16 trials/ecc/cue

% parameters
task{1}.parameter.cue = stimulus.parameter.cue; %0=neutral; 1=exogenous cue
task{1}.parameter.targetEcc = stimulus.gabor.ecc;
task{1}.parameter.targetHemi = [2 1]; % [right left]
task{1}.parameter.targetInterval = stimulus.targetInterval; % [first second]
task{1}.random = 1;

%% Initialize broken trial parameters
stimulus.brokenTrial = [];
stimulus.brokenTrial.orgNumTrials = task{1}.numTrials;
stimulus.brokenTrial.trialIdx = [];
stimulus.noise.noiseEcc = [];
stimulus.noise.noiseTheta = [];

%% Run tasks
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,...
    @screenUpdateCallback,@trialResponseCallback,@startTrialCallback,...
    @endTrialCallback);

% run task
tnum = 1;
while (tnum<=length(task)) && ~myscreen.userHitEsc
    [task, myscreen, tnum] = updateTask(task,myscreen,tnum);
    myscreen = tickScreen(myscreen,task);
end

% display that the next block is loading
mglTextSet([],15,[0 0 0],[],[],[],0);
if myscreen.blocks.n==myscreen.blocks.nBlocks
    message = 'Saving data...';
else
    message = 'Please take a short break before the next block starts.';
end
mglClearScreen; mglFlush; mglClearScreen;
mglTextDraw(message,[0 0]);
mglWaitSecs(1);
mglFlush;
mglWaitSecs(1.5);

endTask(myscreen,task);

%% startTrial
function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

% set stimVal (i.e., contrast)
stimulus.gabor.contrast = stimulus.stair.xCurrent;
stimulus.stimVal(task.trialnum) = stimulus.gabor.contrast;

% keep track of eye position on every frame then average when determining
% if fixation was broken or not. this will hopefully stop the trials to be
% broken on every other trial because the eye tracker keeps thinking my eye
% is jumping around the screen (>1 degree jumps)
stimulus.eye.avgDistFromFix = [];

% store previous zero position. this will be the new zero point of the
% fixation distance calculation
if task.trialnum>1 && stimulus.eye.eyeTrack && ~stimulus.brokenTrial.trialIdx(task.trialnum-1)
    stimulus.eye.prevZeroPos = stimulus.eye.zeroPos;
else
    stimulus.eye.prevZeroPos = [0 0]; % reset to calibration grid center
end
stimulus.eye.trialStart = mglGetSecs; % this is here to stop the initial drift correct from being skipped

%% Initialize index for whether the trial is broken or not
stimulus.brokenTrial.trialIdx(task.trialnum) = 0;

%% Cues
% determine direction of cue (for target)
thetas = [pi 0]; %pi=left; 0=right
targTheta = thetas(task.thistrial.targetHemi);

% in noise interval, target cue is never repeated (use a different
% combination of cue and eccentricity)
cueComb = combvec(thetas,unique(stimulus.gabor.ecc))';
if task.thistrial.targetEcc==0
    cueComb = setdiff(cueComb,[thetas' zeros(length(thetas),1)],'rows');
else
    cueComb = setdiff(cueComb,[targTheta task.thistrial.targetEcc],'rows');
end
cueComb = cueComb(randi(length(cueComb)),:);
noiseTheta = cueComb(1);
noiseEcc = cueComb(2);

% store cued eccentricity and hemifield
stimulus.noise.noiseEcc(task.trialnum) = noiseEcc;
stimulus.noise.noiseTheta(task.trialnum) = noiseTheta;

% determine cue position for both exogenous cue and location cue
[cueX, cueY] = pol2cart([targTheta noiseTheta],[task.thistrial.targetEcc noiseEcc]);
% location cue
stimulus.locCue.x0 = cueX+stimulus.fixation.loc(1);
stimulus.locCue.x1 = stimulus.locCue.x0;
stimulus.locCue.y0 = cueY+stimulus.fixation.loc(2);
stimulus.locCue.y1 = stimulus.locCue.y0;

% Exo
if task.thistrial.cue
    % create the cue
    stimulus.exoCue.x0 = cueX-stimulus.exoCue.rMin; %[target noise]
    stimulus.exoCue.x1 = stimulus.exoCue.x0+stimulus.exoCue.rMax;
    stimulus.exoCue.y0 = cueY+stimulus.exoCue.yDistFromTarg;
    stimulus.exoCue.y1 = stimulus.exoCue.y0;
end


%% startSegment
function [task, myscreen] = startSegmentCallback(task,myscreen)
global stimulus

stimulus.eye.noFixT0 = []; % no-fixation timer
stimulus.eye.fixT0 = []; % fixation timer

switch task.thistrial.thisseg
    case 1
        makeTarget(task)
    case 2
        mglTextSet([],15,[0 0 0],[],[],[],0);
    case 3
        mglTextSet([],20,[0 0 0],[],[],[],1);
    case {4 10}
        if stimulus.eye.eyeTrack
            % take average eye posotion during the fixation period (segment 2)
            % and set as new zero poition (dyanmic fixation updating)
            stimulus.eye.zeroPos = nanmean(stimulus.eye.fixationPos,1);
            stimulus.eye.fixationPos = []; % initialize fixation variable
        end
end

%% screenUpdate
function [task, myscreen] = screenUpdateCallback(task,myscreen)
global stimulus
mglClearScreen

if stimulus.eye.eyeTrack
    % continuously update eye position
    stimulus.eye.currentPos = mglEyelinkGetCurrentEyePos;
    if ~isempty(stimulus.eye.zeroPos)
        dist = calcDistance(stimulus.eye.currentPos,stimulus.eye.zeroPos);
    end
end

switch task.thistrial.thisseg
    case 1
        % texture made, start trial
        if stimulus.target.made
            task = jumpSegment(task);
        end
        % make fixation cross
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
    case 2
        % check if drift correction is necessary
        if stimulus.eye.eyeTrack
            if isempty(stimulus.eye.zeroPos)
                mglTextDraw(['Please fixate on the center of the cross and press ',...
                    'RETURN.'],[0 1]);
                % to get a more precise measurement of the subject's
                % central fixation point, collect eye position throughout
                % this period and when spacebar is pressed, use the mean as
                % the true estimate of center
                stimulus.eye.fixationPos(end+1,:) = stimulus.eye.currentPos;
                
                % if spacebar was pressed
                if mglGetKeys(37) && isempty(stimulus.eye.zeroPos) && ...
                        stimulus.eye.eyeTrack && mglGetSecs(stimulus.eye.trialStart)>=0.1
                    % first check if subject was looking within some range of the center
                    % of the screen, as defined by average fixation position
                    stimulus.eye.fixationPos = nanmean(stimulus.eye.fixationPos,1);
                    if calcDistance(stimulus.eye.fixationPos,stimulus.eye.prevZeroPos)<stimulus.eye.breakTrial
                        stimulus.eye.zeroPos = stimulus.eye.fixationPos;
                    else
                        stimulus.eye.fixationPos = [];
                    end
                end
            else
                % check that subject's eye position is at fixation
                if isempty(stimulus.eye.fixT0)
                    stimulus.eye.fixT0 = mglGetSecs; % initialize counter
                    stimulus.eye.fixationPos = [];
                end
                if dist<=stimulus.eye.breakTrial
                    % update fixation counter
                    fixCounter = mglGetSecs(stimulus.eye.fixT0);
                    
                    % reset non-fixation counter
                    stimulus.eye.noFixT0 = [];
                    
                    if fixCounter>=stimulus.eye.fixTimeReq
                        stimulus.eye.fixT0 = [];
                        task = jumpSegment(task);
                    end
                else
                    % reset counter if subject breaks fixation to
                    % ensure that subject fixates for at least 50 ms
                    stimulus.eye.fixT0 = [];
                    
                    % start counter for non-fixation time
                    if isempty(stimulus.eye.noFixT0)
                        stimulus.eye.noFixT0 = mglGetSecs;
                    end
                    noFixCounter = mglGetSecs(stimulus.eye.noFixT0);
                    % ask subject to drift correct if no fixation occurs in
                    % 3 seconds.
                    if noFixCounter>=stimulus.eye.noFixLimit
                        stimulus.eye.zeroPos = [];
                    end
                end
            end
        end
        % display fixation cross
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
    case {3 9} % fixation period
        if stimulus.eye.eyeTrack
            % collect eye position when fixation is being held
            stimulus.eye.fixationPos(end+1,:) = stimulus.eye.currentPos;
        end
        % display fixation cross
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
    case 4 % first cue
        if task.thistrial.cue
            if task.thistrial.targetInterval==1
                % draw target cue
                mglLines2(stimulus.exoCue.x0(1),stimulus.exoCue.y0(1),stimulus.exoCue.x1(1),...
                    stimulus.exoCue.y1(1),stimulus.exoCue.size,stimulus.exoCue.color);
            elseif task.thistrial.targetInterval==2
                % draw noise cue
                mglLines2(stimulus.exoCue.x0(2),stimulus.exoCue.y0(2),stimulus.exoCue.x1(2),...
                    stimulus.exoCue.y1(2),stimulus.exoCue.size,stimulus.exoCue.color);
            end
        else
            % draw neutral cue
            % upper line
            mglLines2(stimulus.frame.x0,stimulus.frame.y,...
                stimulus.frame.x1,stimulus.frame.y,...
                stimulus.exoCue.size,stimulus.exoCue.color);
            
            % lower line
            mglLines2(stimulus.frame.x0,-stimulus.frame.y,...
                stimulus.frame.x1,-stimulus.frame.y,...
                stimulus.exoCue.size,stimulus.exoCue.color);
        end
        % display fixation cross (because I find my eyes wandering)
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
    case 10 % second cue
        if task.thistrial.cue
            if task.thistrial.targetInterval==2
                % draw target cue
                mglLines2(stimulus.exoCue.x0(1),stimulus.exoCue.y0(1),stimulus.exoCue.x1(1),...
                    stimulus.exoCue.y1(1),stimulus.exoCue.size,stimulus.exoCue.color);
            elseif task.thistrial.targetInterval==1
                % draw noise cue
                mglLines2(stimulus.exoCue.x0(2),stimulus.exoCue.y0(2),stimulus.exoCue.x1(2),...
                    stimulus.exoCue.y1(2),stimulus.exoCue.size,stimulus.exoCue.color);
            end
        else
            % draw neutral cue
            % upper line
            mglLines2(stimulus.frame.x0,stimulus.frame.y,...
                stimulus.frame.x1,stimulus.frame.y,...
                stimulus.exoCue.size,stimulus.exoCue.color);
            
            % lower line
            mglLines2(stimulus.frame.x0,-stimulus.frame.y,...
                stimulus.frame.x1,-stimulus.frame.y,...
                stimulus.exoCue.size,stimulus.exoCue.color);
        end
        % display fixation cross (because I find my eyes wandering)
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
    case {5 7 11 13 15} % short fixation periods
        % display fixation cross
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
    case 6 % interval 1
        if task.thistrial.targetInterval==1
            mglBltTexture(stimulus.target.texture,[0 0]); % target
        else
            mglBltTexture(stimulus.noise.texture,[0 0]); % noise
        end
    case 8 % location cue for interval 1
        if task.thistrial.targetInterval==1 % target
            % draw upper line
            mglLines2(stimulus.locCue.x0(1),...
                stimulus.locCue.y0(1)+stimulus.exoCue.yDistFromTarg,stimulus.locCue.x1(1),...
                stimulus.locCue.y1(1)+stimulus.exoCue.yDistFromTarg+stimulus.locCue.length,...
                stimulus.locCue.size,stimulus.locCue.color);
            % draw lower line
            mglLines2(stimulus.locCue.x0(1),...
                stimulus.locCue.y0(1)-stimulus.exoCue.yDistFromTarg,stimulus.locCue.x1(1),...
                stimulus.locCue.y1(1)-stimulus.exoCue.yDistFromTarg-stimulus.locCue.length,...
                stimulus.locCue.size,stimulus.locCue.color);
        elseif task.thistrial.targetInterval==2 % cue random location
            % draw upper line
            mglLines2(stimulus.locCue.x0(2),...
                stimulus.locCue.y0(2)+stimulus.exoCue.yDistFromTarg,stimulus.locCue.x1(2),...
                stimulus.locCue.y1(2)+stimulus.exoCue.yDistFromTarg+stimulus.locCue.length,...
                stimulus.locCue.size,stimulus.locCue.color);
            % draw lower line
            mglLines2(stimulus.locCue.x0(2),...
                stimulus.locCue.y0(2)-stimulus.exoCue.yDistFromTarg,stimulus.locCue.x1(2),...
                stimulus.locCue.y1(2)-stimulus.exoCue.yDistFromTarg-stimulus.locCue.length,...
                stimulus.locCue.size,stimulus.locCue.color);
        end
        
        % display fixation cross
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
    case 12 % interval 2
        if task.thistrial.targetInterval==2
            mglBltTexture(stimulus.target.texture,[0 0]); % target
        else
            mglBltTexture(stimulus.noise.texture,[0 0]); % noise
        end
    case 14 % location cue for interval 2
        if task.thistrial.targetInterval==2 % target
            % draw upper line
            mglLines2(stimulus.locCue.x0(1),...
                stimulus.locCue.y0(1)+stimulus.exoCue.yDistFromTarg,stimulus.locCue.x1(1),...
                stimulus.locCue.y1(1)+stimulus.exoCue.yDistFromTarg+stimulus.locCue.length,...
                stimulus.locCue.size,stimulus.locCue.color);
            % draw lower line
            mglLines2(stimulus.locCue.x0(1),...
                stimulus.locCue.y0(1)-stimulus.exoCue.yDistFromTarg,stimulus.locCue.x1(1),...
                stimulus.locCue.y1(1)-stimulus.exoCue.yDistFromTarg-stimulus.locCue.length,...
                stimulus.locCue.size,stimulus.locCue.color);
        elseif task.thistrial.targetInterval==1 % cue random location
            % draw upper line
            mglLines2(stimulus.locCue.x0(2),...
                stimulus.locCue.y0(2)+stimulus.exoCue.yDistFromTarg,stimulus.locCue.x1(2),...
                stimulus.locCue.y1(2)+stimulus.exoCue.yDistFromTarg+stimulus.locCue.length,...
                stimulus.locCue.size,stimulus.locCue.color);
            % draw lower line
            mglLines2(stimulus.locCue.x0(2),...
                stimulus.locCue.y0(2)-stimulus.exoCue.yDistFromTarg,stimulus.locCue.x1(2),...
                stimulus.locCue.y1(2)-stimulus.exoCue.yDistFromTarg-stimulus.locCue.length,...
                stimulus.locCue.size,stimulus.locCue.color);
        end
        
        % display fixation cross
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
    case 16
        % make black outline
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size,stimulus.fixation.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size),stimulus.fixation.color);
        
        % make fixation cross green
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            stimulus.fixation.size-stimulus.fixation.edgeOffset,stimulus.exoCue.color);
        mglFillRect(stimulus.fixation.loc(1),stimulus.fixation.loc(2),...
            fliplr(stimulus.fixation.size)-stimulus.fixation.edgeOffset,stimulus.exoCue.color);
end

if stimulus.eye.eyeTrack && ismember(task.thistrial.thisseg,3:length(task.thistrial.seglen)-1)
    % 1) Sometimes the eyetracker just randomly says that I blink,
    % inserting a single nan into my eye position data. So, to avoid
    % breaking trials for inaccurate tracking, at least two nans are
    % required for a blink to be recognized. (I did some testing to confirm
    % that a real blink produces two consecutive frames of nans at 100 Hz.
    
    % break trial if eye position is not at fixation any time between cue
    % presentation and stimulus presentation
    breakTrial = 0;
    blinkThresh = 8;
    if ~ieNotDefined('dist')
        stimulus.eye.avgDistFromFix = [stimulus.eye.avgDistFromFix; dist];
        if sum(stimulus.eye.avgDistFromFix==inf)>=blinkThresh && all(stimulus.eye.avgDistFromFix(end:-1:end-(blinkThresh-1))==inf)% 1
            % if there's a blink, break trial
            breakTrial = 1;
            disp('blink break');
        elseif length(stimulus.eye.avgDistFromFix)>=stimulus.eye.breakTrialTime
            % average across the last n frames (n is specified by
            % breakTrialTime) and use that average distance to determine if
            % fixation was broken
            avgDist = stimulus.eye.avgDistFromFix(end:-1:end-(stimulus.eye.breakTrialTime-1));
            avgDist(avgDist==inf) = nan;
            avgDist = nanmean(avgDist);
            if avgDist>stimulus.eye.breakTrial
                breakTrial = 1;
                disp('dist break');
            end
        end
        
        % break the trial if needed
        if breakTrial
            stimulus.brokenTrial.trialIdx(task.trialnum) = 1;
            task = jumpSegment(task,inf);
        end
    end
end

%% trialResponse
function [task, myscreen] = trialResponseCallback(task,myscreen)
global stimulus

% give auditory feedback if wrong
if find(task.thistrial.buttonState)~=task.thistrial.targetInterval
    mglPlaySound(stimulus.sound.incorrect); % incorrect
    correct = 0;
else
    mglPlaySound(stimulus.sound.correct); % correct
    correct = 1;
end

updateStair = 0;
if task.thistrial.targetEcc<3 && ~task.thistrial.cue
    % only update staircase when the target is within 3 dva.
    updateStair = 1;
end

if updateStair
    stimulus.stair = usePalamedesStaircase(stimulus.stair,correct);
end

% go to next segment
task = jumpSegment(task);

%% endTrialCallback
function [task, myscreen] = endTrialCallback(task,myscreen)
global stimulus

% store segment lengths for each trial to aid eyetracking analysis (so that
% we can isolate segments...etc)
stimulus.eye.seglen(task.trialnum,:) = task.thistrial.seglen;

% store drift-corrected (0,0) point
if stimulus.eye.eyeTrack
    stimulus.eye.driftCorrection(task.trialnum,:) = stimulus.eye.zeroPos;
end

if stimulus.brokenTrial.trialIdx(task.trialnum)
    % play sound for broken fixation
    mglPlaySound('Submarine');
    
    % keep track of the parameters on this trial and then move on to the
    % next trial
    parameters = fieldnames(task.block(task.blocknum).parameter);
    for p = 1:length(parameters)
        if ~isfield(stimulus.brokenTrial,parameters{p})
            stimulus.brokenTrial.(parameters{p}) = [];
        end
        stimulus.brokenTrial.(parameters{p})(end+1) = ...
            task.block(task.blocknum).parameter.(parameters{p})(task.blockTrialnum);
    end
    
    % add a trial to the end of the block
    task.numTrials = task.numTrials+1;
    
    % pause experiment to make sure subject is ready for the next trial
    mglClearScreen; mglFlush;
    mglWaitSecs(1.5);
end

% if we've reached the end of the original number of trials, start
% presenting trials that need to be re-done
if task.trialnum>=stimulus.brokenTrial.orgNumTrials
    parameters = fieldnames(task.block(task.blocknum).parameter);
    
    % if done with the original trials then update the block info to show a
    % division between original and re-done trials
    if task.trialnum==stimulus.brokenTrial.orgNumTrials
        task.block(task.blocknum).trialn = task.blockTrialnum;
        
        % make a new block for the re-done trials
        if isfield(stimulus.brokenTrial,(parameters{1}))
            task.blocknum = task.blocknum+1;
            task.block(task.blocknum).trialn = inf;
            task.blockTrialnum = 0;
        end
    end
    
    % for trials that need to be re-done, continuously update the
    % parameters in that block
    if task.block(task.blocknum).trialn==inf && ~isempty(stimulus.brokenTrial.(parameters{1}))
        for p = 1:length(parameters)
            % replace block structure with parameters of the trial that will be
            % re-done
            if ~isfield(task.block(task.blocknum).parameter,parameters{p})
                task.block(task.blocknum).parameter.(parameters{p}) = [];
            end
            task.block(task.blocknum).parameter.(parameters{p}) = ...
                [task.block(task.blocknum).parameter.(parameters{p}) ...
                stimulus.brokenTrial.(parameters{p})];
            
            % clear that trial, showing that the trial has been re-done
            stimulus.brokenTrial.(parameters{p}) = [];
        end
    end
end

% If n consecutive broken fixations happen, re-do the drift correct
nToReDoDriftCorr = 3;
if task.trialnum>=nToReDoDriftCorr && ...
        sum(stimulus.brokenTrial.trialIdx(end-(nToReDoDriftCorr-1):end))==nToReDoDriftCorr
    stimulus.eye.zeroPos = [];
end

function makeTarget(task)
global stimulus

stimulus.target.made = 0;
if isfield(stimulus,'target.texture')
    mglDeleteTexture(stimulus.target.texture);
    mglDeleteTexture(stimulus.noise.texture);
end

% create texture
temp = make2ndOrderTex(stimulus.carrier.sf,stimulus.gabor.sf,...
    stimulus.gabor.contrast,task.thistrial.targetEcc,...
    task.thistrial.targetHemi,1,'gaborSD',stimulus.gabor.sd,'gaborPhase',...
    stimulus.gabor.phase,'width',stimulus.carrier.width,...
    'height',stimulus.carrier.height);
noiseTemp = make2ndOrderTex(stimulus.carrier.sf,stimulus.gabor.sf,...
    stimulus.gabor.contrast,task.thistrial.targetEcc,...
    task.thistrial.targetHemi,0,'gaborSD',stimulus.gabor.sd,...
    'gaborPhase',stimulus.gabor.phase,'width',stimulus.carrier.width,...
    'height',stimulus.carrier.height);

% format matrix to display faster (as instructed on MGL How Tos)
temp = repmat(uint8(temp),[1 1 3]);
temp = permute(temp,[3 2 1]);
temp(4,:,:) = 255;

noiseTemp = repmat(uint8(noiseTemp),[1 1 3]);
noiseTemp = permute(noiseTemp,[3 2 1]);
noiseTemp(4,:,:) = 255;

temp = mglCreateTexture(temp);
noiseTemp = mglCreateTexture(noiseTemp);
stimulus.target.texture = temp; clear temp;
stimulus.noise.texture = noiseTemp; clear noiseTemp

% set flag to dispaly texture
stimulus.target.made = 1;