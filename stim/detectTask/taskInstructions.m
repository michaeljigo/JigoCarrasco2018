% Allows the subject to step through the frames of the experiment with
% instructions on how to do the task and what to expect.

function taskInstructions

mglClose;
myscreen.displayname = 'gdm';
myscreen.background = [0.5 0.5 0.5];
initScreen(myscreen);
mglTextSet([],15,[0 0 0]);

% start with response keys
mglTextDraw('First, thank you for being a willing participant in this study!',[0 3.5]);
mglTextDraw('For this study you will be using three keys:',[0 1.5]);
mglTextDraw('1 (on number pad): target-present response',[0 0.5]);
mglTextDraw('2 (on number pad): target-absent response',[0 -0.5]);
mglTextDraw('RETURN: advance through instructions.',[0 -1.5]);
mglTextDraw('Press RETURN to continue',[0 -3.5]);
mglFlush

while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end
mglClearScreen;

% show target and noise stimuli
orgDir = pwd;
cd('./textures/');

% show the gabor
x = make2ndOrderTex(2,0.25,1,1,1,2,'dispFig=0','newCarrier=1');
x = mglCreateTexture(x);
mglBltTexture(x,[0 0]);
mglTextDraw('Your target will be this pattern:',[0 6]); mglFlush;
while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end

% show carrier
mglClearScreen;
x = make2ndOrderTex(2,0.25,1,1,1,0,'dispFig=0','newCarrier=1');
x = mglCreateTexture(x);
mglBltTexture(x,[0 0]);
mglTextDraw('Embedded in this noisy image.',[0 6]); mglFlush
while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end

% show modulated carrier (THRICE)
commandwindow
mglClearScreen; start = mglGetSecs;
mglTextDraw('Which will make the target look like this (press RETURN to see different targets):',[0 6]);
x1 = make2ndOrderTex(2,0.25,1,1,1,1,'dispFig=0','newCarrier=1');
x2 = make2ndOrderTex(2,0.25,1,1,1,1,'dispFig=0','newCarrier=1');
x3 = make2ndOrderTex(2,0.25,1,1,1,1,'dispFig=0','newCarrier=1');
x4 = make2ndOrderTex(2,0.25,1,1,1,1,'dispFig=0','newCarrier=1');
x5 = make2ndOrderTex(2,0.25,1,1,1,1,'dispFig=0','newCarrier=1');
x1 = mglCreateTexture(x1);
x2 = mglCreateTexture(x2);
x3 = mglCreateTexture(x3);
x4 = mglCreateTexture(x4);
x5 = mglCreateTexture(x5);
mglBltTexture(x1,[0 0]); mglFlush; pause

mglClearScreen;
mglTextDraw('Which will make the target look like this (press RETURN to see different targets):',[0 6]);
mglBltTexture(x2,[0 0]); mglFlush; pause

mglClearScreen;
mglTextDraw('Which will make the target look like this (press RETURN to see different targets):',[0 6]);
mglBltTexture(x3,[0 0]); mglFlush; pause

mglClearScreen;
mglTextDraw('Which will make the target look like this (press RETURN to see different targets):',[0 6]);
mglBltTexture(x4,[0 0]); mglFlush; pause

mglClearScreen;
mglTextDraw('Which will make the target look like this (press RETURN to see different targets):',[0 6]);
mglBltTexture(x5,[0 0]); mglFlush; pause


% show target at different eccentricities
% left hemi
left = make2ndOrderTex(2,0.25,1,7.2,1,1,'dispFig=0','newCarrier=1');
% right hemi
right = make2ndOrderTex(2,0.25,1,7.2,2,1,'dispFig=0','newCarrier=1');
mglClearScreen;
mglTextDraw('The location of the target will vary on each trial between these two extremes (left):',[0 6]);
left = mglCreateTexture(left);
mglBltTexture(left,[0 0]); mglFlush; pause; mglClearScreen;

mglTextDraw('The location of the target will vary on each trial between these two extremes (right):',[0 6]);
right = mglCreateTexture(right);
mglBltTexture(right,[0 0]); mglFlush; pause; mglClearScreen;

% show possible target locations
mglTextDraw('At the beginning of each block, you will be shown the possible target locations.',[0 8]);
mglTextDraw('Importantly, each location is tagged with a number. You should remember these number-location associations.',[0 7]);
mglTextDraw('There will be one of two sets of target locations in each block. I call this the "large spacing" set',[0 6]);
stimulus.gabor.ecc = [0 2.4 4.8 7.2];
x = make2ndOrderTex(2,0.25,1,7.2,2,0,'dispFig=0','newCarrier=1');
x = mglCreateTexture(x);
mglBltTexture(x,[0 0]);

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
mglFlush; pause; mglClearScreen;

mglTextDraw('And this is the "small spacing" set',[0 6]);
stimulus.gabor.ecc = [0 1.2 3.6 6];
x = make2ndOrderTex(2,0.25,1,7.2,2,0,'dispFig=0','newCarrier=1');
x = mglCreateTexture(x);
mglBltTexture(x,[0 0]);

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
mglFlush;

cd(orgDir); pause; mglClearScreen;



% Now go frame by frame in experiment
start = mglGetSecs; doneWithScreen = 0;
stimulus.fixation.width = 0.4;
stimulus.fixation.size = 1.2;
stimulus.fixation.color = [0 0 0];
stimulus.fixation.loc = [-0.02 0.02];
mglTextDraw('In each trial, you will have to fixate on a central cross.',[0 7]);
mglTextDraw('If you move your eyes away from this cross, the trial will end and the following sound will play...',[0 6]);
mglTextDraw('Press RETURN to play the sound.',[0 5]);
mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
    stimulus.fixation.color,stimulus.fixation.loc);
mglFlush;

while ~doneWithScreen
    while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end
    if mglGetKeys(37) && mglGetSecs(start)>0.5
        mglPlaySound('Submarine');
        mglWaitSecs(1); mglClearScreen;
        doneWithScreen=1;
    end
end

% Cue
stimulus.carrier.width = 30.5; %30.5
stimulus.carrier.height = 10; %10
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

mglTextDraw('When you fixate, a cue will appear. The cue can either be a frame',[0 6]);
% draw frame
mglPolygon(stimulus.frame.outerX,stimulus.frame.outerY,...
    stimulus.fixation.color);
mglPolygon(stimulus.frame.innerX,stimulus.frame.innerY,...
    myscreen.background);

% display fixation cross
mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
    stimulus.fixation.color,stimulus.fixation.loc);
mglFlush; pause; mglClearScreen;

% Peripheral cue
stimulus.cue.theta = 0; task.thistrial.targetEcc = 2.4;
mglTextDraw('A peripheral line segment indicating the target location',[0 6]);
stimulus.gabor.sd = 0.6;
stimulus.exoCue.rMin = 0.5; % starting point of line
stimulus.exoCue.rMax = 1; % extent of line
stimulus.exoCue.yDistFromTarg = 2*stimulus.gabor.sd; % this ensures that the bar will be above 95% of the gabor
% create the cue
[stimulus.exoCue.x0, stimulus.exoCue.y0] = pol2cart(stimulus.cue.theta,task.thistrial.targetEcc);
stimulus.exoCue.x0 = stimulus.exoCue.x0-stimulus.exoCue.rMin;
stimulus.exoCue.x1 = stimulus.exoCue.x0+stimulus.exoCue.rMax;
stimulus.exoCue.y0 = stimulus.exoCue.y0+stimulus.exoCue.yDistFromTarg;
stimulus.exoCue.y1 = stimulus.exoCue.y0;
mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
    stimulus.fixation.color,stimulus.fixation.loc);
mglLines2(stimulus.exoCue.x0,stimulus.exoCue.y0,stimulus.exoCue.x1,...
    stimulus.exoCue.y1,stimulus.fixation.size,stimulus.fixation.color);
mglFlush; pause; mglClearScreen;

% Central cue
mglTextDraw('Or a central line segment with a number telling you how far from fixation the target will appear',[0 7]);
mglTextDraw('You should use this information to help you do the task.',[0 6]);
stimulus.respCue.dispNum = 1;
stimulus.respCue.rMin = 0.25;
stimulus.respCue.rMax = 0.55;
stimulus.respCue.xNum = -0.04; stimulus.respCue.yNum = 0;
if task.thistrial.targetEcc>0
    [stimulus.respCue.x0, stimulus.respCue.y0] = pol2cart(stimulus.cue.theta,...
        stimulus.respCue.rMin);
    [stimulus.respCue.x1, stimulus.respCue.y1] = pol2cart(stimulus.cue.theta,...
        stimulus.respCue.rMax);
end
% display number at origin
mglTextDraw(num2str(stimulus.respCue.dispNum),...
    [stimulus.respCue.xNum,stimulus.respCue.yNum]);
if stimulus.respCue.dispNum>0
    mglLines2(stimulus.respCue.x0,stimulus.respCue.y0,...
        stimulus.respCue.x1,stimulus.respCue.y1,...
        stimulus.fixation.size,[0 0 0]);
end
mglFlush; pause; mglClearScreen;

% Target presentation
mglTextDraw('Then the noisy stimulus will be briefly displayed',[0 7]);
mglTextDraw('The target will be embedded in this stimulus on half of the trials.',[0 6]);
cd('./textures');
right = make2ndOrderTex(2,0.25,1,2.4,2,1,'dispFig=0','newCarrier=1');
cd(orgDir);
right = mglCreateTexture(right);
mglBltTexture(right,[0 0]);
mglFixationCross(stimulus.fixation.width,stimulus.fixation.size,...
    stimulus.fixation.color,stimulus.fixation.loc);
mglFlush; pause; mglClearScreen;

% Response cue
mglTextDraw('Following the stimulus, a number and line segment will appear.',[0 7]);
mglTextDraw('The number and line segment asks the following question: Did the target appear at indicated location?',[0 6]);
mglTextDraw('This is when you can press 1 or 2 to report whether or not the target was present. Feel free to blink when you respond.',[0 5]);
mglTextDraw(num2str(stimulus.respCue.dispNum),[stimulus.respCue.xNum,stimulus.respCue.yNum]);
if stimulus.respCue.dispNum>0
    mglLines2(stimulus.respCue.x0,stimulus.respCue.y0,stimulus.respCue.x1,...
        stimulus.respCue.y1,stimulus.fixation.size,[0 0 0]);
end
mglFlush; pause; mglClearScreen;

% Auditory feedback
stimulus.sound.incorrect = mglInstallSound('/users/purplab/Desktop/jigo/experiments/asot/stim/incorrect.wav');
stimulus.sound.correct = mglInstallSound('/users/purplab/Desktop/jigo/experiments/asot/stim/correct.wav');
mglTextDraw('When you respond correctly, you will hear this sound (please press RETURN)',[0 0]); mglFlush;
doneWithScreen = 0;
while ~doneWithScreen
    while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end
    if mglGetKeys(37) && mglGetSecs(start)>0.5
        mglPlaySound(stimulus.sound.correct);
        mglWaitSecs(1); mglClearScreen;
        doneWithScreen=1;
    end
end

mglTextDraw('When you are incorrect, you will hear this sound (please press RETURN)',[0 0]); mglFlush;
doneWithScreen = 0;
while ~doneWithScreen
    while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end
    if mglGetKeys(37) && mglGetSecs(start)>0.5
        mglPlaySound(stimulus.sound.incorrect);
        mglWaitSecs(1); mglClearScreen;
        doneWithScreen=1;
    end
end

mglTextDraw('Now that you know the structure of the task, you will do some practice to get used to the weird stimulus.',[0 0]); mglFlush;
mglWaitSecs(5); mglClose;