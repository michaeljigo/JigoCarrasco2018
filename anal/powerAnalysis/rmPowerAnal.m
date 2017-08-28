% perform a power analysis for 2x2 repeated measures designs

% adapted from
% https://cognitivedatascientist.com/2015/12/14/power-simulation-in-r-the-repeated-measures-anova-5/

% each column of MU is the first independent variable (within-subject
% factor) and each row of MU is the second within-subject factor

function [avgSim, powerRes] = rmPowerAnal(data,nSubj,nRep)

mu = data.avg;
sigma = data.ste;
rho = corr(data.avg(:,1),data.avg(:,2));

% make the covariance matrix
sigma = repmat(sigma,2,2);
covarMat = sigma.*sigma*rho;
covarMat(logical(eye(2))) = unique(sigma).^2;

% simulate data
for s = 1:length(nSubj)
    disppercent(-inf/nRep,sprintf('Simulating data...'));
    for rep = 1:nRep
        % initialize matrix that will hold simulated data
        subj = nan(nSubj(s),size(mu,2),size(mu,1));
        % loop through each level of the second within-subject factor
        for l = 1:size(mu,1)
            levelMu = mu(l,:);
            % generate data for each subject
            subj(:,:,l) = mvnrnd(levelMu,covarMat,nSubj(s));
            
            if all(mu<1)
                % constrain performance between 0 and 1
                subj(subj>1) = 1;
            end
            subj(subj<0) = 0;
        end
        
        % do repeated-measures ANOVA on simulated data
        temp = buildrmaov2Matrix(subj);
        temp = rmaov2(temp,'printOut=0');
        
        % store p-values for main effects and interactions
        iv1(rep) = temp.iv1{3};
        iv2(rep) = temp.iv2{3};
        iv12(rep) = temp.iv12{3};
        disppercent(rep/nRep);
        
        avgSim(rep,:,:) = squeeze(mean(subj,1));
    end
    % store power results
    powerRes.power.iv1(s) = sum(iv1<0.05)/nRep;
    powerRes.power.iv2(s) = sum(iv2<0.05)/nRep;
    powerRes.power.iv12(s) = sum(iv12<0.05)/nRep;
    disppercent(inf,sprintf('Simulation %i complete',s));
end
% store the number of subjects used in simulation
powerRes.n = nSubj;