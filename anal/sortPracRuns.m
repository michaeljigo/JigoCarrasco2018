% sort practice data into a folder named 'prac' for each subject

function sortPracRuns(expName)

switch expName
    case 'carL_gabL_exo'
        subj = {'AJF' 'AF' 'CB' 'CP' 'MJ' 'MR'};
    case 'carL_gabL_endo'
        subj = {'AF' 'AJF' 'CB' 'CP' 'MJ' 'MR'};
    otherwise
        fprintf('Cute experiment name.\n');
        return
end

for s = 1:length(subj)
    % first go into the folder and make the prac folder if non-existent, if
    % it does exist, remove it to start fresh and avoid duplicates
    pracDir = ['../data/',subj{s},'/',expName,'/prac/'];
    if ~exist(pracDir,'dir')
        mkdir(pracDir);
    else
        rmdir(pracDir,'s');
        mkdir(pracDir);
    end
    
    % now move threshold files into the prac folder
    threshDir = ['../data/',subj{s},'/',expName,'/thresh/'];
    try
        copyfile([threshDir,'*stim*.mat'],pracDir);
    catch
    end
    
    % now search through doNotInclude for old files that were deemed as
    % practice and put them into the practice folder in the main directory
    try
        dniDir = '../data/doNotInclude/';
        dni = dir([dniDir,'/',subj{s},'*']);
        dniDir = [dniDir,dni.name,'/',expName,'/prac/'];
        copyfile([dniDir,'*stim*'],pracDir);
    catch
        fprintf('%s has no extra practice files \n',subj{s});
    end
end