% runs endogenous cueing texture segmentation task with second-order
% textures

function runHighContrast(myscreen)
global stimulus % initStimulus in main driver function

orgStimulus = stimulus;

stimulus.eye.blinkBreak = 0;
stimulus.eye.distBreak = 0;

%% Set parameters for stimuli
% cue
stimulus.respCue.rMin = 0.25;
stimulus.respCue.rMax = 0.55;

% fixation
stimulus.fixation.width = 0.4;
stimulus.fixation.size = 1.2;
stimulus.fixation.color = [0 0 0];
stimulus.fixation.loc = [-0.02 0.02]; % for some reason the fixation cross isn't perfectly centered, this fixes it

%% Set up task parameters
eyeContingent = stimulus.eye.eyeContingent; % screen to make sure subject is fixating
fixation = 0.5;
cue = stimulus.cue.dur;
isi = stimulus.soa;
stim = 0.1; % 10 frames (larger than maximum duration from Yeshurun & Carrasco, 2000)
respIsi = 0.2;
respCue = inf;
minIti = 0.1;
maxIti = 0.6; % the real iti = eyeContingent+iti;

task{1}.segmin = [eyeContingent fixation cue isi stim respIsi respCue minIti];
task{1}.segmax = [eyeContingent fixation cue isi stim respIsi respCue maxIti];
task{1}.getResponse = [0 0 0 0 0 0 1 0];
task{1}.waitForBacktick = 0;
task{1}.numTrials = 28; %arbritrarily set for now

% parameters
task{1}.parameter.targetEcc = stimulus.gabor.ecc;
task{1}.parameter.targetHemi = [2 1]; % [right left]
task{1}.parameter.targetPresent = [0 1]; % [absent present]
task{1}.random = 1;
task{1}.randVars.uniform.carrierFlip = 1:size(stimulus.tex.textures,1);

%% Initialize broken trial parameters
stimulus.brokenTrial = [];
stimulus.brokenTrial.orgNumTrials = task{1}.numTrials;
stimulus.brokenTrial.trialIdx = [];

%% Initialize stimulus variables
stimulus.stimVal = [];

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
% before ending the experiment, delete the textures from the stimulus
% structure because it takes up quite a bit of space
stimulus = rmfield(stimulus,'tex');

% display that the next block is loading
mglClearScreen; mglFlush; mglClearScreen;
mglTextDraw('Please wait for the next block.',[0 0]);
mglWaitSecs(1);
mglFlush;
mglWaitSecs(3);

stimulus = orgStimulus;
% myscreen.saveData = 0;
% endTask(myscreen,task);

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

%% Index pre-made texture
if isfield(stimulus.tex,'thisTexture')
    mglDeleteTexture(stimulus.tex.thisTexture);
end
if task.thistrial.targetPresent
    % index texture to be displayed during the trial
    h = find(stimulus.tex.texParams.hemifield==task.thistrial.targetHemi);
    e = find(stimulus.tex.texParams.eccentricity==task.thistrial.targetEcc);
    c = find(stimulus.tex.texParams.contrast==1);
else
    % index the last element in the first dimension, this will be the noise
    % patch
    h = size(stimulus.tex.textures,2);
    e = 1;
    c = 1;
end
% get carrier flip state
f = task.thistrial.carrierFlip;
stimulus.tex.thisTexture = mglCreateTexture(stimulus.tex.textures{f,h,e,c});

%% Response cue
switch task.thistrial.targetHemi
    case 2 % rightward along horizontal meridian
        stimulus.cue.theta = 0;
    case 1 % leftward along horizontal meridian
        stimulus.cue.theta = pi;
end
% only display arrow if cue is not 0
if task.thistrial.targetEcc>0
    [stimulus.cue.x0, stimulus.cue.y0] = pol2cart(stimulus.cue.theta,...
        stimulus.respCue.rMin);
    [stimulus.cue.x1, stimulus.cue.y1] = pol2cart(stimulus.cue.theta,...
        stimulus.respCue.rMax);
end

% display number at origin
stimulus.cue.xNum = -0.04; stimulus.cue.yNum = 0;
stimulus.cue.dispNum = find(unique(stimulus.gabor.ecc)==task.thistrial.targetEcc);
stimulus.cue.dispNum = stimulus.cue.dispNum-1;

%% Neutral cue
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

%% startSegment
function [task, myscreen] = startSegmentCallback(task,myscreen)
global stimulus

stimulus.eye.noFixT0 = []; % no-fixation timer
stimulus.eye.fixT0 = []; % fixation timer

if task.thistrial.thisseg==3 && stimulus.eye.eyeTrack
    stimulus.eye.zeroPos = nanmean(stimulus.eye.fixationPos,1);
    stimulus.eye.fixationPos = []; % initialize fixation variable
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
                    % ensure that subject fixates for at least 100 ms
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
    case 2 % fixation period
        if stimulus.eye.eyeTrack
            % collect eye position when fixation is being held
            stimulus.eye.fixationPos(end+1,:) = stimulus.eye.currentPos;
        end
        % display fixation cross
        mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
            stimulus.fixation.color,stimulus.fixation.loc);
    case {4 6 8}
        % display fixation cross
        mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
            stimulus.fixation.color,stimulus.fixation.loc);
    case 3 % cue
        % draw frame
        mglPolygon(stimulus.frame.outerX,stimulus.frame.outerY,...
            stimulus.fixation.color);
        mglPolygon(stimulus.frame.innerX,stimulus.frame.innerY,...
            myscreen.background);
        
        % display fixation cross
        mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
            stimulus.fixation.color,stimulus.fixation.loc);
    case 5 % texture
        mglBltTexture(stimulus.tex.thisTexture,[0 0]);
    case 7 % response cue
        % display number at origin
        mglTextDraw(num2str(stimulus.cue.dispNum),[stimulus.cue.xNum,stimulus.cue.yNum]);
        if stimulus.cue.dispNum>0
            mglLines2(stimulus.cue.x0,stimulus.cue.y0,stimulus.cue.x1,...
                stimulus.cue.y1,stimulus.fixation.size*0.9,[0 0 0]);
        end
end


if stimulus.eye.eyeTrack && ismember(task.thistrial.thisseg,2:6)
    % 1) Sometimes the eyetracker just randomly says that I blink,
    % inserting a single nan into my eye position data. So, to avoid
    % breaking trials for inaccurate tracking, at least two nans are
    % required for a blink to be recognized. (I did some testing to confirm
    % that a real blink produces two consecutive frames of nans at 100 Hz.
    
    % break trial if eye position is not at fixation any time between cue
    % presentation and stimulus presentation
    breakTrial = 0;
    blinkThresh = 8; % blinkThresh nans in a row are required for a blink to register
    if ~ieNotDefined('dist')
        stimulus.eye.avgDistFromFix = [stimulus.eye.avgDistFromFix; dist];
        if sum(stimulus.eye.avgDistFromFix==inf)>=blinkThresh && all(stimulus.eye.avgDistFromFix(end:-1:end-(blinkThresh-1))==inf)% 1
            % if there's a blink, break trial
            breakTrial = 1; disp('blink break');
            stimulus.eye.blinkBreak = stimulus.eye.blinkBreak+1;
        elseif length(stimulus.eye.avgDistFromFix)>=stimulus.eye.breakTrialTime
            % average across the last n frames (n is specified by
            % breakTrialTime) and use that average distance to determine if
            % fixation was broken
            avgDist = stimulus.eye.avgDistFromFix(end:-1:end-(stimulus.eye.breakTrialTime-1));
            avgDist(avgDist==inf) = nan;
            avgDist = nanmean(avgDist);
            if avgDist>stimulus.eye.breakTrial
                breakTrial = 1; disp('dist break');
                stimulus.eye.distBreak = stimulus.eye.distBreak+1;
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
if task.thistrial.buttonState(1)~=task.thistrial.targetPresent
    sound(sin(0.5*pi*(0:pi/32:75)));
else
    sound(sin(1.5*pi*(0:pi/32:75)));
end

% update staircase with hits and misses
stimulus.stimVal(task.trialnum) = stimulus.stair.xCurrent;
updateStair = 0;

if task.thistrial.targetPresent && task.thistrial.targetEcc<6
    % if a target is presented, only update the staircase when the
    % subject responds correctly or incorrectly to eccentricities
    % within 6 dva
    if task.thistrial.buttonState(1)
        hit = 1; % hit
    else
        hit = 0; % miss
    end
    updateStair = 1;
end

if updateStair
    stimulus.stair = usePalamedesStaircase(stimulus.stair,hit);
    % if the staircase is stuck on the maximum or minimum value for 10 trials
    % in a row, reset to beginning
    stuckTrials = 10;
    if length(stimulus.stair.x)>stuckTrials
        stuckTrials = stuckTrials-1;
        maxStim = max(stimulus.stair.priorAlphaRange);
        minStim = min(stimulus.stair.priorAlphaRange);
        if all(abs(stimulus.stair.x(end-stuckTrials:end)-maxStim)<1*10^-10) || ...
                all(abs(stimulus.stair.x(end-stuckTrials:end)-minStim)<1*10^-10)
            % reset posterior and xCurrent
            stimulus.stair.pdf = stimulus.stair.prior;
            [stimulus.stair.xCurrent, ~, ~] = PAL_AMRF_pdfDescriptives(stimulus.stair.pdf,...
                stimulus.stair.priorAlphaRange);
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