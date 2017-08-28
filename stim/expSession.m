% INPUTS:
% expSession(subj,expPrefix,eyeTrack,initializeSubj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use to run the experimental session. ensures that endo and exo
% session order is counterbalanced across subjects

% expName = car[LH]_gab[LH]

function expSession(subj,expPrefix,eyeTrack,initializeSubj)

% initialize experiment variables
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
    thisEccOrder = expOrder(subjIdx).toDo.ecc{currentIdx}(1);
end

if ~initializeSubj
    fprintf('====== %s cue (eccSet #%i) ======\n',thisCue,thisEccOrder);
    main_texDetect_withRespCue(subj,0,thisEccOrder,8,[expPrefix,'_',thisCue],eyeTrack);
    
    % after the experiment is done update the subject list
    expOrder(subjIdx).done.ecc{currentIdx}(end+1) = expOrder(subjIdx).toDo.ecc{currentIdx}(1);
    expOrder(subjIdx).toDo.ecc{currentIdx}(1) = [];
    
    if isempty(expOrder(subjIdx).toDo.ecc{currentIdx})
        expOrder(subjIdx).done.cue{currentIdx} = expOrder(subjIdx).master.cue{currentIdx};
        expOrder(subjIdx).toDo.cue{currentIdx} = [];
    end
    save(filename,'subjList','expOrder');
else
    % just initialize subject without running the tasks
    save(filename,'subjList','expOrder');
    expOrder(subjIdx).toDo.cue
end