function funWrap(expName)

switch expName
    case 'carL_gabL_exo'
        subj = {'AJF' 'AF' 'CB' 'CP' 'MJ' 'MR' 'KN' 'HL'};
    case 'carL_gabL_endo'
        subj = {'AF' 'AJF' 'CB' 'CP' 'MJ' 'MR' 'KN' 'HL'};
    otherwise
        fprintf('Cute experiment name.\n');
        return
end

doFit.y = 1;
doFit.fittype = 'polynomial';
doFit.grpCoeff = [];

for s = 1:length(subj)
    anal_cue(subj{s},expName,1,0,doFit); close all;
%             eyePosition(subj{s},expName);
%     anal_dPrime(subj{s},expName,1,doFit); close all
% bootstrapPeaks_dPrime(subj{s},expName); pause; close all
end