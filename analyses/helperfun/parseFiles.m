% Parses all the stimfiles in a directory into sets that does (but not
% always, depending on the value of nSets) contain at least one of each
% condition

% Input:
% dataDir = path to directory

% nSets = number of sets to create.
% extractFields = the fields within the exp structure (output from
% getTaskParameters) from which data will be extracted

% conditionStr = the variable that holds the condition labels. For example,
% if stimulus.attn is a value of 0 or 1 corresponding to exogenous or
% endogenous cues, conditionVar would be stimulus.attn
% The input is a string.

% extractFields = cell array of strings specifying the fields in the
% stimulus or in the task structure that you want to extract data from

% Output:
% parsedFiles is a structure array containing the data from the files
% stimFiles contains the names of the files that were loaded

function [parsedFiles, stimFiles] = parseFiles(dataDir,nSets,extractFields,varargin)

% set defaults of task and phase numbers
getArgs(varargin,{'taskNum=[]','phaseNum=[]','conditionStr=nan'});

% set default condition string, such that all extracted variables will be
% under the same condition label, i.e., 1.
if isnan(conditionStr)
    conditionStr = '0';
end

% get the files
if ~strcmp(dataDir(end),'/')
    dataDir = [dataDir,'/'];
end
files = dir([dataDir,'*stim*.mat']);

if isempty(files)
    error('no files to load');
end

allFiles = []; parsedFiles = [];
%% First...
% Go through the files and separate files into their respective cueing
% conditions
for i = 1:length(files)
    load([dataDir,files(i).name]);
    % make condition field to store the number of conditions being analyzed
    evalc(['conditionVar =', conditionStr]);
    if ~isfield(allFiles,'condition')
        allFiles.condition = conditionVar;
        totalFiles(1) = 0;
    end
    
    % check if the condition label already exists in the condition field,
    % if not add it
    if ~ismember(allFiles.condition,conditionVar)
        allFiles.condition(end+1) = conditionVar;
        totalFiles(end+1) = 0;
    end
    
    % specify the value that will index the condition, the value (i.e., 1,
    % 2, or 3) will be mapped onto the position of the condition label in
    % allFiles.condition, such that cond1 holds the data for the condition
    % in the first element of allFiles.condition
    condIdx = find(ismember(allFiles.condition,conditionVar));
    
    % now put data file names into the corresponding fields named by
    % condIdx
    if ~isfield(allFiles,['cond',num2str(condIdx)])
        allFiles.(['cond',num2str(condIdx)]) = [];
    end
    allFiles.(['cond',num2str(condIdx)]){end+1} = files(i).name;
    
    % count the total number of files for each condition
    totalFiles(condIdx) = totalFiles(condIdx)+1;
end

% now separate files into sets
if nSets==inf % an infinite input tells the program to keep all stimfiles separate
    nSets = totalFiles;
end
filesPerSet = max(ceil(totalFiles/nSets));

% give a warning message if the inputted number of sets cannot be possibly
% made
if filesPerSet*nSets>max(totalFiles)
    fprintf('Not enough files to make %i sets. Making %i instead. \n',nSets,...
        ceil(max(totalFiles)/filesPerSet));
end

maxFiles = find(max(totalFiles)==totalFiles); % use largest number of files to determine when to stop parsing files
maxFiles = maxFiles(1); % arbritrarily used the first index
totalSets = ceil(totalFiles(maxFiles)/filesPerSet);
setN = 0; % initialize counter for number of sets

% sort the conditions for purely superficial reasons
[~,sortIdx] = sort(allFiles.condition);
temp.condition = allFiles.condition(sortIdx);
for i = 1:length(sortIdx)
    temp.(['cond',num2str(i)]) = allFiles.(['cond',num2str(sortIdx(i))]);
end
allFiles = temp;
clearvars -except allFiles dataDir extractFields conditionStr taskNum ...
    phaseNum totalSets filesPerSet

%% SECOND (Initialize the parsedFiles structure)
% Go through each condition, load a file, get the name of each condition we
% want to extract
parsedFiles = []; holdingFields = [];
for c = 1:length(allFiles.condition)
    % load file
    load([dataDir,allFiles.(['cond',num2str(c)]){1}]);
    exp = getTaskParameters(myscreen,task);
    
    % get data from the specified task and phase
    if isempty(taskNum)
        taskIdx = '(1)';
    else
        taskIdx = sprintf('{%i}',taskNum);
    end
    if isempty(phaseNum)
        exp = eval(['exp',taskIdx]);
    else
        phaseIdx = sprintf('(%i)',phaseNum);
        exp = eval(['exp',taskIdx,phaseIdx]);
    end
    
    % get variables
    for v = 1:length(extractFields)
        switch extractFields{v}
            case {'randVars' 'parameter'}
                f = fieldnames(exp.(extractFields{v}));
                for iF = 1:length(f)
                    if isfield(parsedFiles,f{iF})
                        continue
                    else
                        parsedFiles.(f{iF}) = [];
                        % store, in a separate varaible, the fields that hold
                        % the desired variables
                        holdingFields{end+1} = ['exp.',extractFields{v},...
                            '.',f{iF}];
                    end
                end
            otherwise
                %if ~isempty(parsedFiles) && any(cellfun(@(x) ...
                        %~isempty(strfind(extractFields{v},x)),fieldnames(parsedFiles)))
                if ~isempty(parsedFiles) && ismember(extractFields{v},fieldnames(parsedFiles))
                    continue
                else
                    if isfield(exp,extractFields{v})
                        parsedFiles.(extractFields{v}) = [];
                        holdingFields{end+1} = ['exp.',extractFields{v}];
                    else
                        % assumed to be stimulus
                        try
                        temp = eval(['stimulus.',extractFields{v}]);
                        clear temp
                        % set the fieldName to be used in parsedFiles to be
                        % the characters after the last period
                        fieldName = strfind(extractFields{v},'.');
                        if isempty(fieldName)
                            parsedFiles.(extractFields{v}) = [];
                        else
                            fieldName = extractFields{v}(fieldName(end)+1:end);
                            parsedFiles.(fieldName) = [];
                        end
                        holdingFields{end+1} = ['stimulus.',extractFields{v}];
                        catch
                            warning([extractFields{v}, ' is not a real variable.']);
                        end
                    end
                end
        end
    end
end
% add in condition field
parsedFiles.condition = [];
% make structure array reflect number of sets
parsedFiles = repmat(parsedFiles,1,totalSets);
% make structure that will hold stim files that go into each set
stimFiles.name = [];
stimFiles = repmat(stimFiles,1,totalSets);

clearvars -except allFiles dataDir extractFields conditionStr taskNum ...
    phaseNum totalSets parsedFiles filesPerSet holdingFields stimFiles

%% THIRD (Add in the variables into parsedFiles structure)
parsedFields = fieldnames(parsedFiles);
parsedFields = setdiff(parsedFields,'condition');

% Create counter to keep track of parsing
nTotal = totalSets*length(allFiles.condition)*filesPerSet*length(parsedFields);
nCount = 0;
%disppercent(-inf/nTotal,sprintf('Parsing files...'));

for s = 1:totalSets
    % determine number of files that CAN be put into each set from each
    % condition
    startFileNum = filesPerSet*(s-1)+1;
    endFileNum = startFileNum+filesPerSet-1;
    for c = 1:length(allFiles.condition)
        for f = startFileNum:endFileNum
            % if the condition has less than the number of files that CAN
            % be put into a set, skip the condition
            if f>length(allFiles.(['cond',num2str(c)]))
                fprintf('Condition #%i has less files than expected. Ignored.\n',c);
                break
            end
            stimFiles(s).name{end+1} = allFiles.(['cond',num2str(c)]){f};
            % load data and extract the intended variables
            load([dataDir,allFiles.(['cond',num2str(c)]){f}]);
            exp = getTaskParameters(myscreen,task);
            
            % get data from the specified task and phase
            if isempty(taskNum)
                taskIdx = '(1)';
            else
                taskIdx = sprintf('{%i}',taskNum);
            end
            if isempty(phaseNum)
                exp = eval(['exp',taskIdx]);
            else
                phaseIdx = sprintf('(%i)',phaseNum);
                exp = eval(['exp',taskIdx,phaseIdx]);
            end
            
            % loop through the initialized fields and add in the variables
            for h = 1:length(parsedFields)
                % determine where to get the data from
                data = cellfun(@(x) ~isempty(strfind(x,parsedFields{h})),...
                    holdingFields);
                % check if the desired field exists in the current data
                % file
                try 
                    data = eval(holdingFields{data});
                    
                    % verify that the same length of vector is inserted into
                    % the parsedFiles structure...
                    % to do this, use a field with exp (i.e., the output from
                    % getTaskParameters) to determine the number of elements
                    % that need to be added
                    %nElements = find(cellfun(@(x) ~isempty(strfind(x,'exp')),...
                        %holdingFields));
                    %nElements = min(nElements);
                    %nElements = length(eval(holdingFields{nElements}));
                    %data = data(1:nElements);
                    parsedFiles(s).(parsedFields{h}) = ...
                        horzcat(parsedFiles(s).(parsedFields{h}),data);
                catch
                    warning([allFiles.(['cond',num2str(c)]){f}, ' does not',...
                        ' contain the field: ',holdingFields{data},'. ',...
                        'Skipping file.']);
                    continue
                end
                
                % Update counter
                nCount = nCount+1;
                %disppercent(nCount/nTotal);
            end
        end
        % After each condition, pad structure to compensate for non-existent
        % variable(s) in the current condition or for unequal variables
        % (parameters) across conditions, which would be bad experimental
        % design
        desiredLen = max(structfun(@(x) length(x),parsedFiles(s)));
        parsedFiles(s) = structfun(@(x) validateVarLen(x,desiredLen), ...
            parsedFiles(s),'UniformOutput',false);
        
        % Add in condition value
        parsedFiles(s).condition(isnan(parsedFiles(s).condition)) = ...
            allFiles.condition(c);
    end
end
%disppercent(inf,sprintf('Files have been parsed, and it only'));

function out = validateVarLen(in,desiredLen)
if length(in)~=desiredLen
    out = nan(1,desiredLen-length(in));
    out = [in out];
else
    out = in;
end
