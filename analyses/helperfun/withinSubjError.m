function err = withinSubjError(data,morey)
   % Usage:    withinErr = withinSubjErr_morey(data)
   % Purpose:  Program for sizing error bars in within-subjects experimental designs. Assumes
   %           data is set up with the first dimension as a subject and each following 
   %           dimension as a fixed, within-subject factor. 
   % By:       Michael Jigo, 08/15/19
   %
   %           References: Morey 2008, "Confidence intervals from normalized data: 
   %                       A correction to Cousineau (2005)"
   %
   %                       Cousineau 2005, "Confidence intervals in within-subject designs:
   %                       A simpler solution to Loftus and Masson's method"
   %
   % Input:    data    --  nxm matrix where n is the # of subjects and m is the # of conditions 
   %                       data can be of any dimensionality. the first dimension just needs to 
   %                       correspond to individual subjects
   %           morey   --  flag for Morey's correction. 0=no correction; 1=Morey correction (default) 
   % 
   % Output:   err     --  within-subject standard error of the mean for each condition

if ~exist('morey','var')
   morey = 1;
end

% compute necessary variables...
dim = size(data);
nSubj = dim(1); nCond = dim(2:end); 
subjMeans = nanmean(data(:,:),2);
grandMean = nanmean(data(:));

% normalize data w.r.t subject and group mean
data = data-repmat(subjMeans,[1 nCond]); % remove subject mean
data = data+grandMean; % add grand mean

% compute new sample error of the mean
err = squeeze(nanstd(data,[],1)./sqrt(nSubj));

% apply Morey's correction
if morey
   M = prod(nCond);
   M = sqrt(M/(M-1));
   err = squeeze(err*M);
end
return;

