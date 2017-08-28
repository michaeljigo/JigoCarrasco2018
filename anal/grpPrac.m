% determine how many practice trials subjects performed for each task

function grpPrac(exclude)

if ieNotDefined('exclude')
    exclude = 0;
end
nBlocks = [];
exps = {'carL_gabL_endo' 'carL_gabL_exo'};
for e = 1:length(exps)
    groupDir = ['../data/group/',exps{e},'/prac/'];
    
    % where subject results are stored
    subj = dir ([groupDir,'*.mat']);
    
    % get subject initials to make sure correct subjects are being used in
    % group analysis
    subjNames = arrayfun(@(x) x.name(1:strfind(x.name,'.')-1), subj,...
        'UniformOutput',false);
    % display subject names
    disp('%%%%% Subjects %%%%%');
    cellfun(@(x) fprintf('%s \n',x),subjNames);
    if exclude
        try
            subjExc = eval(input(sprintf('Enter the subjects you want excluded (use cell array of strings): '),'s'));
            subjExc = ismember(subjNames,subjExc);
            subj(subjExc) = [];
            subjNames(subjExc) = [];
        catch
            fprintf('Not removing any subjects.\n');
        end
    end
    
    % loop through subjects and load data
    for s = 1:length(subj)
        load([groupDir,subj(s).name]);
        
        nBlocks(end+1) = (prac.neutral+prac.cue)/56;
    end
    
end
% just print the mean and standard deviation
fprintf('Avg.: %.1f +- %.1f blocks \n',mean(nBlocks),std(nBlocks));