% INPUTS:
% expSession(subj,expPrefix,eyeTrack,workspaceTex,initializeSubj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use to run the experimental session. ensures that endo and exo
% session order is counterbalanced across subjects

% expName = car[LH]_gab[LH]

function expSession(subj,expPrefix,eyeTrack,workspaceTex,initializeSubj)

% initialize experiment variables
if ieNotDefined('workspaceTex')
    workspaceTex = [];
end
if ieNotDefined('eyeTrack')
    eyeTrack = 1;
end
if ieNotDefined('initializeSubj')
    initializeSubj = 0;
end

% first check that there is a subject list for this experiment saved in
% the current directory
filename = ['./.subjList_',expPrefix,'.mat'];
if exist(filename,'file')
    % load the file
    load(filename);
else
    % create a new subject list
    subjList = {subj};
    eccOrder = [1 2]; % set the possible eccentricity order
    eccOrder = eccOrder(randperm(length(eccOrder))); % randomize eccentricity set
    
    expOrder(1).master.cue = {'exo' 'endo'};
    expOrder(1).master.ecc = {eccOrder fliplr(eccOrder)}; % counterbalance within subjects the eccentricity set
    expOrder(1).toDo.cue = expOrder(1).master.cue;
    expOrder(1).toDo.ecc = expOrder(1).master.ecc;
    expOrder(1).done.cue = cell(1,length(expOrder(1).master.cue));
    expOrder(1).done.ecc = expOrder(1).done.cue;
end

% check if the subject is already present in the subjList
[~,subjIdx] = ismember(subj,subjList);

if all(subjIdx==0)
    % if subject is not on the list, add them and set the order based on
    % the order of the previous subject to ensure counterbalancing
    subjList{end+1} = subj;
    subjIdx = length(subjList);
    
    eccOrder = [1 2];
    eccOrder = eccOrder(randperm(length(eccOrder)));
    
    expOrder(end+1).master.cue = fliplr(expOrder(end).master.cue);
    expOrder(end).master.ecc = {eccOrder fliplr(eccOrder)};
    expOrder(end).toDo.cue = expOrder(end).master.cue;
    expOrder(end).toDo.ecc = expOrder(end).master.ecc;
    expOrder(end).done.cue = cell(1,length(expOrder(end).master.cue));
    expOrder(end).done.ecc = expOrder(end).done.cue;
end

% check what cue the subject will use in this session
if all(cellfun('isempty',expOrder(subjIdx).toDo.cue))
    % if both cue conditions are done
    fprintf('%s is done with both experiments.\n',subj);
    return
else
    % get current cue
    currentIdx = min(find(~cellfun('isempty',expOrder(subjIdx).toDo.cue)));
    thisCue = expOrder(subjIdx).toDo.cue{currentIdx};
    % get the order of eccentricity set
    thisEccOrder = expOrder(subjIdx).toDo.ecc{currentIdx};
end

if ~initializeSubj
    for e = 1:length(thisEccOrder)
        fprintf('====== %s cue (eccSet #%i) ======\n',thisCue,thisEccOrder(e));
        myscreen = main_texDetect(subj,0,thisEccOrder(e),6,[expPrefix,'_',thisCue],eyeTrack,workspaceTex); 
        % let subject rest and then have a quick refresher of the target
        % locations
        if e==1
            eccentricityRefresher(thisEccOrder,thisCue,myscreen,workspaceTex);
        end
    end
    
    % after the experiment is done update the subject list
    expOrder(subjIdx).done.ecc{currentIdx} = expOrder(subjIdx).toDo.ecc{currentIdx};
    expOrder(subjIdx).toDo.ecc{currentIdx} = [];
    
    if isempty(expOrder(subjIdx).toDo.ecc{currentIdx})
        expOrder(subjIdx).done.cue{currentIdx} = expOrder(subjIdx).master.cue{currentIdx};
        expOrder(subjIdx).toDo.cue{currentIdx} = [];
    end
    save(filename,'subjList','expOrder');
    
    % display thank you message
    mglClearScreen;
    mglTextDraw('Thank you for participating!',[0 0]);
    mglTextDraw('Please let the experimenter know you are done.',[0 -1]);
    mglFlush;
    mglWaitSecs(5);
    mglClose;
else
    % just initialize subject without running the tasks
    save(filename,'subjList','expOrder');
    expOrder(subjIdx).toDo.cue
end


function eccentricityRefresher(thisEccOrder,thisCue,myscreen,workspaceTex)
global stimulus

mglClearScreen;
mglTextDraw('You are done with the first half of this session!',[0 1]);
mglTextDraw('If you want to take a quick break, please do so. When you are ready to continue press RETURN.',[0 0]);
mglFlush;
% wait for button press to move on
mglWaitSecs(0.5);
while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end
mglClearScreen; mglFlush; mglClearScreen;

if thisEccOrder(2)==1
    setName = 'large spacing';
    eccSet = [0 repmat([2.4 4.8 7.2],1,2)];
elseif thisEccOrder(2)==2
    setName = 'small spacing';
    eccSet = [0 repmat([1.2 3.6 6],1,2)];
end

mglTextDraw(['For the next half of the session, you will be using the ',...
    upper(setName), ' set shown below.'],[0 6]);

% draw carrier noise that will be used during each trial
% create texture to be used during trial
tex = workspaceTex.textures{1,3,1,1};
tex = mglCreateTexture(tex);
mglBltTexture(tex,[0 0]);

% draw the 0 first
mglTextDraw('0',[0 0]);
ecc = unique(eccSet);
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
mglTextDraw(['To get you used to this new set, ',...
    'you will do 14 trials with these new target locations.'],[0 -5.5]);
mglTextDraw('Press RETURN to start those trials.',[0 -6.5]);
mglFlush;
% wait for button press to move on
mglWaitSecs(0.5);
while ~any(mglGetKeys(37)); mglWaitSecs(0.1); end
mglClearScreen; mglFlush; mglClearScreen;

% run high contrast task with new eccentricity set
stimulus.gabor.contrast = 1;
stimulus.gabor.ecc = eccSet;
stimulus.targetTilt = randi(2);
stimulus.numTrials = 14;
stimulus.parameter.cue = 1; % only do trials with the cue
stimulus.tex.textures = workspaceTex.textures;
stimulus.tex.texParams = workspaceTex.texParams;
myscreen.saveData = 0;
stimulus.eye.eyeTrack = 0;
stimulus.eye.eyeContingent = 0.01;
if strcmp(thisCue,'exo')
    runExoTexDetect(myscreen);
elseif strcmp(thisCue,'endo')
    runEndoTexDetect(myscreen);
end
clear global stimulus myscreen