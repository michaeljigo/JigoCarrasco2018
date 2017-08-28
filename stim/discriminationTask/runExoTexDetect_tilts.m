% runs exogenous cueing texture segmentation task with second-order
% textures

function runExoTexDetect_tilts(myscreen)
global stimulus % initStimulus in main driver function

%% Set parameters for stimuli

stimulus.stimVal = [];
%1) this ensures that the bar will be 0.3 degrees above one standard
% deviation of gabor

% fixation
stimulus.fixation.width = 0.4;
stimulus.fixation.size = 1.2;
stimulus.fixation.color = [0 0 0];
stimulus.fixation.loc = [-0.02 0.02]; % for some reason the fixation cross isn't perfectly centered, this fixes it

% response cue
stimulus.respCue.rMin = 0.2;
stimulus.respCue.rMax = 0.45;

% exogenous cue
stimulus.exoCue.rMin = 0.35; %0.35 starting point of line
stimulus.exoCue.rMax = stimulus.exoCue.rMin*2; % extent of line
stimulus.exoCue.yDistFromTarg = stimulus.gabor.sd+(stimulus.gabor.sd/2); % 1
% this variable lets me run subjects through a quick refresher session
% between eccentricity sets
if ~isfield(stimulus,'targetTilt')
    stimulus.targetTilt = [1 2]; %[1=CCW; 2=CW]
    stimulus.numTrials = 112;
end

%% Set up task parameters
eyeContingent = stimulus.eye.eyeContingent; % screen to make sure subject is fixating
fixation = 0.5;
cue = 0.04; %0.04
cueISI = 0.05; %0.05
stim = 0.1; % 10 frames (larger than maximum duration from Yeshurun & Carrasco, 2000)
resp = 0.3;
respISI = 0.2;
makeTex = inf; % make texture at the beginning of a trial, and its duration will be included in the iti

task{1}.seglen = [makeTex eyeContingent fixation cue cueISI stim respISI resp];
task{1}.getResponse = [zeros(1,length(task{1}.seglen)-1) 1];
task{1}.waitForBacktick = 0;
task{1}.numTrials = stimulus.numTrials; % 16 trials/ecc/cue

% parameters
task{1}.parameter.cue = stimulus.parameter.cue; %0=neutral; 1=exogenous cue
task{1}.parameter.targetEcc = stimulus.gabor.ecc;
task{1}.parameter.targetHemi = [2 1]; % [right left]
task{1}.parameter.targetTilt = stimulus.targetTilt; % [vertical tilted]
task{1}.parameter.tiltDir = [1 2]; % [CCW CW]
task{1}.random = 1;

%% Initialize broken trial parameters
stimulus.brokenTrial = [];
stimulus.brokenTrial.orgNumTrials = task{1}.numTrials;
stimulus.brokenTrial.trialIdx = [];

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
mglClearScreen; mglFlush; mglClearScreen;
if myscreen.blocks.n==myscreen.blocks.nBlocks
    if strcmp(myscreen.subj,'JP')
        message = 'I fixed the message for you Jake. You happy now princess?';
    else
        message = 'Saving data...';
    end
else
    message = 'Please take a short break before the next block starts.';
end
mglTextDraw(message,[0 0]);
mglWaitSecs(1);
mglFlush;
mglWaitSecs(3);

endTask(myscreen,task);

%% startTrial
function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

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

% Exo
if task.thistrial.cue
    % create the cue
    [cueX, cueY] = pol2cart(targTheta,task.thistrial.targetEcc);
    stimulus.exoCue.x0 = cueX-stimulus.exoCue.rMin; %[target noise]
    stimulus.exoCue.x1 = stimulus.exoCue.x0+stimulus.exoCue.rMax;
    stimulus.exoCue.y0 = cueY+stimulus.exoCue.yDistFromTarg;
    stimulus.exoCue.y1 = stimulus.exoCue.y0;
else
    % Neutral cue
    % prep frame
    x = [stimulus.carrier.width+0.2 stimulus.carrier.width]/2;
    y = [stimulus.carrier.height+0.2 stimulus.carrier.height]/2;
    xOut = x(1); xIn = x(2); yOut = y(1); yIn = y(2);
    frameXScale = [-1 -1 1 1 -1];
    frameYScale = [1 -1 -1 1 1];
    % inner frame
    stimulus.frame.innerX = repmat(xIn,1,length(frameXScale)).*frameXScale;
    stimulus.frame.innerY = repmat(yIn,1,length(frameYScale)).*frameYScale;
    % outer frame
    stimulus.frame.outerX = repmat(xOut,1,length(frameXScale)).*frameXScale;
    stimulus.frame.outerY = repmat(yOut,1,length(frameYScale)).*frameYScale;
end

%% Response cue
% only display arrow if cue is not 0
if task.thistrial.targetEcc>0
    [stimulus.respCue.x0, stimulus.respCue.y0] = pol2cart(targTheta,...
        stimulus.respCue.rMin);
    [stimulus.respCue.x1, stimulus.respCue.y1] = pol2cart(targTheta,...
        stimulus.respCue.rMax);
end
% display number at origin
stimulus.respCue.xNum = -0.04; stimulus.respCue.yNum = 0;
stimulus.respCue.dispNum = find(unique(task.parameter.targetEcc)==task.thistrial.targetEcc);%find(unique(stimulus.gabor.ecc)==task.thistrial.targetEcc);
stimulus.respCue.dispNum = stimulus.respCue.dispNum(1)-1;

%% startSegment
function [task, myscreen] = startSegmentCallback(task,myscreen)
global stimulus

stimulus.eye.noFixT0 = []; % no-fixation timer
stimulus.eye.fixT0 = []; % fixation timer

switch task.thistrial.thisseg
    case 1
        makeTarget(task)
    case 4
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
        mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
            stimulus.fixation.color,stimulus.fixation.loc);
    case 3 % fixation period
        if stimulus.eye.eyeTrack
            % collect eye position when fixation is being held
            stimulus.eye.fixationPos(end+1,:) = stimulus.eye.currentPos;
        end
        % display fixation cross
        mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
            stimulus.fixation.color,stimulus.fixation.loc);
    case 4 % cue
        if task.thistrial.cue
            % draw cue
            mglLines2(stimulus.exoCue.x0,stimulus.exoCue.y0,stimulus.exoCue.x1,...
                stimulus.exoCue.y1,stimulus.fixation.size,stimulus.fixation.color);
        else
            % draw frame
            mglPolygon(stimulus.frame.outerX,stimulus.frame.outerY,...
                stimulus.fixation.color);
            mglPolygon(stimulus.frame.innerX,stimulus.frame.innerY,...
                myscreen.background);
        end
        % display fixation cross (because I find my eyes wandering)
        mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
            stimulus.fixation.color,stimulus.fixation.loc);
    case {5 7}
        % display fixation cross
        mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
            stimulus.fixation.color,stimulus.fixation.loc);
    case 6 % texture
        mglBltTexture(stimulus.target.texture,[0 0]);
    case 8 % response cue
        mglTextDraw(num2str(stimulus.respCue.dispNum),[stimulus.respCue.xNum,stimulus.respCue.yNum]);
        if stimulus.respCue.dispNum>0
            mglLines2(stimulus.respCue.x0,stimulus.respCue.y0,stimulus.respCue.x1,...
                stimulus.respCue.y1,stimulus.fixation.size,stimulus.fixation.color);
        end
end

if stimulus.eye.eyeTrack && ismember(task.thistrial.thisseg,2:length(task.thistrial.seglen)-1)
    % 1) Sometimes the eyetracker just randomly says that I blink,
    % inserting a single nan into my eye position data. So, to avoid
    % breaking trials for inaccurate tracking, at least two nans are
    % required for a blink to be recognized. (I did some testing to confirm
    % that a real blink produces two consecutive frames of nans at 100 Hz.
    
    % break trial if eye position is not at fixation any time between cue
    % presentation and stimulus presentation
    breakTrial = 0;
    blinkThresh = 6;
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
if find(task.thistrial.buttonState)~=task.thistrial.targetTilt
    mglPlaySound(stimulus.sound.incorrect); % incorrect
    correct = 0;
else
    mglPlaySound(stimulus.sound.correct); % correct
    correct = 1;
end

% store staircase stimulus values
% if task.thistrial.targetTilt==1
%     stimulus.stimVal(task.trialnum) = stimulus.stair.ccw.xCurrent;
% elseif task.thistrial.targetTilt==2
%     stimulus.stimVal(task.trialnum) = stimulus.stair.cw.xCurrent;
% end
if task.thistrial.targetTilt==1
    stimulus.stimVal(task.trialnum) = 0;
elseif task.thistrial.targetTilt==2
    if task.thistrial.tiltDir==1
        stimulus.stimVal(task.trialnum) = stimulus.stair.ccw.xCurrent;
    elseif task.thistrial.tiltDir==2
        stimulus.stimVal(task.trialnum) = stimulus.stair.cw.xCurrent;
    end
end

% reset staircase if stuck
if length(stimulus.stair.ccw.x)>=5 && ...
        all(stimulus.stair.ccw.x(end:-1:end-4)==min(stimulus.stair.ccw.priorAlphaRange))
    stimulus.stair.ccw.pdf = stimulus.stair.ccw.prior;
elseif length(stimulus.stair.cw.x)>=5 && ...
        all(stimulus.stair.cw.x(end:-1:end-4)==min(stimulus.stair.cw.priorAlphaRange))
    stimulus.stair.cw.pdf = stimulus.stair.cw.prior;
end
    
if task.thistrial.targetEcc<3 && ~task.thistrial.cue
    %     switch task.thistrial.targetTilt
    %         case 1
    %             % update ccw staircase
    %             stimulus.stair.ccw = usePalamedesStaircase(stimulus.stair.ccw,correct);
    %         case 2
    %             % update cw staircase
    %             stimulus.stair.cw = usePalamedesStaircase(stimulus.stair.cw,correct);
    %     end
    if task.thistrial.targetTilt==2
        if task.thistrial.tiltDir==1
            % update ccw staircase
            stimulus.stair.ccw = usePalamedesStaircase(stimulus.stair.ccw,correct);
        elseif task.thistrial.tiltDir==2
            % update cw staircase
            stimulus.stair.cw = usePalamedesStaircase(stimulus.stair.cw,correct);
        end
    end
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

if stimulus.brokenTrial.trialIdx(task.trialnum)|| ~task.thistrial.gotResponse
    if stimulus.brokenTrial.trialIdx(task.trialnum)
        mglPlaySound('Submarine');
        displayMessage = 'Please fixate';
    else
        mglPlaySound(stimulus.sound.incorrect);
        % set the trial as broken
        stimulus.brokenTrial.trialIdx(task.trialnum) = 1;
        displayMessage = 'Too slow';
    end
    
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
    mglClearScreen; 
    mglTextDraw(displayMessage,[0 0]);
    mglFlush;
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
end

% set gabor angle
% switch task.thistrial.targetTilt
%     case 1
%         % CCW
%         gaborAngle = 90 + stimulus.stair.ccw.xCurrent;
%     case 2
%         % CW
%         gaborAngle = 90 - stimulus.stair.cw.xCurrent;
% end

switch task.thistrial.targetTilt
    case 1
        % vertical
        gaborAngle = 90;
    case 2
        % tilted CCW or CW
        if task.thistrial.tiltDir==1
            % CCW
            gaborAngle = 90 + stimulus.stair.ccw.xCurrent;
        elseif task.thistrial.tiltDir==2
            % CW
            gaborAngle = 90 - stimulus.stair.cw.xCurrent;
        end
end

% create texture
temp = make2ndOrderTex(stimulus.carrier.sf,stimulus.gabor.sf,...
    stimulus.gabor.contrast,task.thistrial.targetEcc,...
    task.thistrial.targetHemi,1,'gaborAngle',gaborAngle,...
    'gaborSD',stimulus.gabor.sd,'gaborPhase',0);
% store texture on each trial
if task.trialnum==1
    stimulus.target.all = nan(task.numTrials,size(temp,1),size(temp,2));
end
stimulus.target.all(task.trialnum,:,:) = temp;
% format matrix to display faster (as instructed on MGL How Tos)
temp = repmat(uint8(temp),[1 1 3]);
temp = permute(temp,[3 2 1]);
temp(4,:,:) = 255;
temp = mglCreateTexture(temp);
stimulus.target.texture = temp; clear temp;

% set flag to dispaly texture
stimulus.target.made = 1;