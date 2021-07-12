% EXAMPLE:
% main_texDetect(subj,thresh,eccSet,nSets,expName,eyeTrack)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepares the screen and some stimulus parameters for
% attention texture segmentation task on second-order textures

% inputs:
% subj: subject initials
% thresh:
% 1 = thresholding session
% 2 = quick threshold check during main experiment days

% eccSet: set of eccentricities to display
% 1 = [0 2.4 4.8 7.2] dva
% 2 = [0 1.2 3.6 6] dva

% nSets: number of sets of blocks
% THRESHOLDING: 1 set = 2 blocks (1 block/eccentricity set)
% MAIN: 1 set = 1 block

% expName: name of experiment and folder in which data will be saved
% car[LH]_gab[LH]_[exo endo] where only one character/word within the brackets is used
% for each experiment

% eyeTrack: use eyetracking?
% 1 = use eyetracker (default)
% 2 = don't use eyetracker

% workspaceTex: if textures are preloaded into the MATLAB environment,
% input the variable name containing the textures here.

function main_texDetect_withRespCue(subj,thresh,eccSet,nSets,expName,eyeTrack)
clear global stimulus

global stimulus
addpath('./textures/');

%% Setup experiment and stimulus parameters

% check experiment name for accuracy
carLH = expName(strfind(expName,'car')+3);
gabLH = expName(strfind(expName,'gab')+3);
exoOrEndo = expName(max(strfind(expName,'_'))+1:end);
if ~ismember(carLH,{'L' 'H'}) || ~ismember(gabLH,{'L' 'H'}) || ...
        ~ismember(exoOrEndo,{'endo' 'exo'})
    error('expName should be car[LH]_gab[LH]_[exo endo]');
end


% carrier noise
if strcmp(carLH,'L')
    stimulus.carrier.sf = 2;
elseif strcmp(carLH,'H')
    stimulus.carrier.sf = 4;
end
stimulus.carrier.width = 30.5; %30.5
stimulus.carrier.height = 10; %10 in Yeshurun & Carrasco, 2000


% gabor
if strcmp(gabLH,'L')
    stimulus.gabor.sf = 0.25;
    stimulus.gabor.phase = 0;
    stimulus.gabor.sd = 0.8;
elseif strcmp(gabLH,'H')
    stimulus.gabor.sf = 1;
    stimulus.gabor.phase = 0;
    stimulus.gabor.sd = 0.6;
end


% load sounds
stimulus.sound.incorrect = mglInstallSound('/users/purplab/Desktop/jigo/experiments/asot/stim/incorrect.wav');
stimulus.sound.correct = mglInstallSound('/users/purplab/Desktop/jigo/experiments/asot/stim/correct.wav');


% eyetracking
if ieNotDefined('eyeTrack')
    stimulus.eye.eyeTrack = 1; % use eyetracker
else
    stimulus.eye.eyeTrack = eyeTrack;
end
if stimulus.eye.eyeTrack
    load('./eyelinkParams.mat');
    myscreen.eyelinkParams = eyelinkParams;
end


% eccentricity settings based on eccSet
if sum(ismember(eccSet,[1 2]))
    % doubled the peripheral eccentricities to have the same number of
    % trials as the foveal location across both eccSets
    eccSet1 = [0 repmat([2.4 4.8 7.2],1,2)];
    eccSet2 = [0 repmat([1.2 3.6 6],1,2)];
else
    error('eccSet should be 1 or 2');
end

% Staircase parameters
stairParams.whichStair = 1; % using best PEST
stairParams.alphaRange = 0.1:0.02:1;
stairParams.fitBeta = 2;
stairParams.fitLambda = 0.01;
stairParams.fitGamma = 0.5;
stairParams.threshPerformance = 0.75;
stairParams.PF = 'arbWeibull';

% if it is thresholding, the number of blocks should always be a multiple
% of 2 because I want to test each set of eccentricities
if thresh   
    % set the SOAs
    switch exoOrEndo
        case 'exo'
            stimulus.soa = 0.05;
            stimulus.cue.dur = 0.04;
        case 'endo'
            stimulus.soa = 0.35;
            stimulus.cue.dur = 0.15;
    end
    switch thresh
        case 1 % staircase session
            nBlocks = nSets*2;
            dataDir = ['../data/',subj,'/',expName,'/thresh/'];
            stimulus.trialnum = 28; %84, should take ~3 minutes/block
        case 2 % threshold check
            nBlocks = nSets*2;
            dataDir = ['../data/',subj,'/',expName,'/threshCheck/'];
            stimulus.trialnum = 56; %84
            % quick threshold check before main experiment
            % load threshold and set prior for staircase
            load([dataDir(1:end-6),'/',subj,'_thresh.mat']);
            % prior will be gaussian, centered on threshold from thresholidng
            % session and have a standard deviation such that 99% of the
            % distribution (i.e., 3 standard deviations) will be within 10
            % steps (each with the magnitude of the minimum step size) of the
            % mean.
            priorStd = min(abs(diff(stairParams.alphaRange)))*40;
            priorStd = priorStd/3;
            prior = normpdf(stairParams.alphaRange,finalThresh,priorStd);
            prior = prior/sum(prior);
            stairParams.lastPosterior = prior;
    end
else % main experiment
    nBlocks = nSets;
    % set the eccentricities to use during the experimental eccSets
    evalc(['stimulus.gabor.ecc = eccSet',num2str(eccSet),';']);
    dataDir = ['../data/',subj,'/',expName,'/cue/'];
    
    % check if subject has already done a block
    prevRun = dir([dataDir,'*.mat']);
    if ~isempty(prevRun)
        % load the posterior as the new starting point
        data = load([dataDir,prevRun(end).name]);
        stairParams.lastPosterior = data.stimulus.stair.pdf;
        stimulus.stair = usePalamedesStaircase(stairParams);
    else
        % load posterior from staircase runs
        loadDir = strfind(dataDir,'/cue/');
        loadDir = dataDir(1:loadDir);
        load([loadDir,'/thresh/',subj,'_thresh.mat']);
        stairParams.lastPosterior = posterior;
        stimulus.stair = usePalamedesStaircase(stairParams);
    end
    stimulus.parameter.cue = [0 1];
end

%% Setup myscreen
myscreen.keyboard.nums = [84 85]; %[1 2] on number pad;
myscreen.displayName = 'gdm';
myscreen.saveData = -2;
myscreen.autoCloseScreen = 0;
myscreen.datadir = dataDir;
if ~exist(myscreen.datadir,'dir')
    mkdir(myscreen.datadir);
end
myscreen.background = [0.5 0.5 0.5];
% initialize the screen
myscreen = initScreen(myscreen);
% initialize the stimulus
myscreen = initStimulus('stimulus',myscreen);

% set text size
mglTextSet([],15,[0 0 0]);

%% Run task
myscreen.blocks.nBlocks = nBlocks;
for n = 1:nBlocks
    myscreen.blocks.n = n;
    % Initialize eyetracker at the beginning of each block
    if stimulus.eye.eyeTrack
        % 1) Trial will break if eye distance exceeds this number (in visual
        % angle)
        
        % 2) If eye position is outside fixation for this period of time,
        % the trial will break.
        
        % 3) Set to infinite to make sure subjects' eye position is truly at
        % fixation before the trial begins.
        
        % 4) Time required for fixation before trial to be considered a real
        % fixation (in seconds)
        
        % 5) Initialize eyetracking drift correct so that subject will drift
        % correct at the beginning of every block
        
        % 6) Limit for how long no-fixation will be tolerated before asking
        % the subject to drift correct again.
        
        % 7) Holds the the zero-position of the eye during trials that were
        % not broken. Initializes at 0,0 for the first drift correction
        % calculation.
        if n>1
            myscreen = eyeCalibDisp(myscreen,sprintf('Block %i/%i, please press RETURN',n,nBlocks));
        else
            myscreen = eyeCalibDisp(myscreen,'If you are ready to start, press RETURN');
        end
        stimulus.eye.breakTrial = 1.5; %1
        stimulus.eye.breakTrialTime = round(0.05/myscreen.frametime); %2
        stimulus.eye.eyeContingent = inf; %3
        stimulus.eye.fixTimeReq = 0.1; %4
        stimulus.eye.fixationPos = [];
        stimulus.eye.zeroPos = []; %5
        stimulus.eye.noFixLimit = 2; %6
        stimulus.eye.prevZeroPos = [0 0]; %7
    else
        stimulus.eye.eyeContingent = 0.01;
    end
    
    % run the experiments
    if thresh
        askIfReady(thresh,exoOrEndo)
        
        % Each half of blocks will have one eccentricity set and have a
        % continuous staircase. The staircase will reset for the second
        % half of blocks.
        if iseven(n)
            % let staircase continue from last block
            stairParams.lastPosterior = stimulus.stair.pdf;
            % increase proportion of 0 and 2.4 eccentricities
            stimulus.gabor.ecc = [0 0 0 2.4 2.4 4.8 7.2];
        elseif isodd(n)
            % initialize staircase
            if thresh==1
                stairParams.lastPosterior = normpdf(stairParams.alphaRange,0.6,0.2);
            end
            stimulus.gabor.ecc = [0 0 0 1.2 1.2 3.6 6];
        end
        stimulus.stair = usePalamedesStaircase(stairParams);
        
        % display target locations
        dispCuedLoc(stimulus)
        
        % run task
        runContrastThresh_withRespCue(myscreen);
    else
        % run the experiment
        askIfReady(thresh,exoOrEndo)
        % display target locations
        dispCuedLoc(stimulus)
        
        if strcmp(exoOrEndo,'endo')
            runEndoTexDetect_withRespCue(myscreen)
        elseif strcmp(exoOrEndo,'exo')
            runExoTexDetect_withRespCue(myscreen)
        end
        
        % continue the staircase from the previous block
        stairParams.lastPosterior = stimulus.stair.pdf;
        stimulus.stair = usePalamedesStaircase(stairParams);
    end
end

% if thresh
    % display thank you message
    mglTextDraw('Thank you for participating!',[0 0]);
    mglTextDraw('Please let the experimenter know you are done.',[0 -1]);
    mglFlush;
%     mglWaitSecs(3);
%     mglClose;
% end

function askIfReady(isThresh,whichCue)

if isThresh>0
    message = 'Training session';
else
    switch whichCue
        case 'exo'
            message = 'Peripheral cue';
        case 'endo'
            message = 'Central cue';
    end
end

mglTextDraw(message,[0 2.5]);
mglTextDraw('1st interval: Press 1',[0 0.5]);
mglTextDraw('2nd interval: Press 2',[0 -0.5]);
mglTextDraw('Press RETURN to start the next block',[0 -3.5]);
mglFlush;

mglWaitSecs(0.5);
while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end
return

function dispCuedLoc(stimulus)
%% Display the cues on the texture
mglTextSet([],20,[0 0 0],[],[],[],0);

% draw carrier noise that will be used during each trial
% create texture to be used during trial
tex = make2ndOrderTex(stimulus.carrier.sf,stimulus.gabor.sf,0,0,1,0,...
    'width',stimulus.carrier.width,'height',stimulus.carrier.height);
tex = mglCreateTexture(tex);
mglBltTexture(tex,[0 0]);

% draw the 0 first
mglTextDraw('0',[0 0]);
ecc = unique(stimulus.gabor.ecc);
ecc = ecc(2:end);

% make numbers into cell of strings
numbers = num2cell(1:length(ecc));
numbers = cellfun(@(x) num2str(x),numbers,'UniformOutput',false);

% now draw all the other eccentricities in all the hemifields
theta = [0 pi];
ecc = repmat(ecc,size(theta,2),1);
theta = repmat(theta',1,size(ecc,2));
[x, y] = pol2cart(theta,ecc);

for e = 1:size(ecc,2)
    for t = 1:size(theta,1)
        mglTextDraw(numbers{e},[x(t,e) y(t,e)]);
    end
end

% put instructions
mglTextSet([],15,[0 0 0],[],[],[],0);
mglTextDraw('Target locations',[0 7]);
mglTextDraw('The target will appear in these locations.',...
    [0,6]);
mglTextDraw('Please press RETURN',[0 -6]);
mglFlush

% wait for button press to move on
mglWaitSecs(0.5);
while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end
return
