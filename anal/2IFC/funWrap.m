function funWrap(expName)

switch expName
    case 'carL_gabL_exo'
        subj = {'MJ' 'AF' 'LR'};
    case 'carL_gabL_endo'
        subj = {'MJ' 'AF' 'XW' 'YJZ' 'LR'};
    otherwise
        fprintf('Cute experiment name.\n');
        return
end

for s = 1:length(subj)
    anal_cue(subj{s},expName); close all;
end