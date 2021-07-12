% Purpose:  Recreate Figure 2. 
%           Results for the Endogenous attention experiment
%
% By:       Michael Jigo
% Edited:   07.01.21
%
% Output:   stats    structure containing results of statistical analyses

function stats = figure2
addpath(genpath('./helperfun'));

   %% Get subject-level results
      % isolate endo data
         [data params]  = pack_experiment_results;
         idx            = find(ismember(params.attntype,'endo')); 
         data           = data(idx);


   %% Compute group-average and error
      % d-prime
         group.avg.dPrime  = squeeze(nanmean(data.dPrime,1));
         group.err.dPrime  = withinSubjError(data.dPrime,0);

      % reaction time
         group.avg.rt      = exp(squeeze(nanmean(log(data.rt),1)));
         group.err.rt      = withinSubjError(data.rt);

      % peak eccentricity
         group.avg.peakEcc = squeeze(nanmean(data.model.peakEcc,1));
         group.err.peakEcc = quantile(squeeze(nanmean(data.boot.peakEcc,1)),[0.025 1-0.025],2); % 95% confidence interval


   %% Perform stats
      % d-prime
         stats.dPrime      = simple_mixed_anova(data.dPrime,[],{'cue' 'ecc'});

      % criterion
         stats.criterion   = simple_mixed_anova(data.criterion,[],{'cue' 'ecc'});
         
      % peak eccentricity
         [~,stats.peakEcc] = ttest(data.model.peakEcc(:,1),data.model.peakEcc(:,2));


   %% Plot
      colors = [0 0 0; 202 0 32]./255;
      figure('name','Endogenous attention','position',[617 165 215 380]);
      
      % discriminability
         subplot(2,1,1)
         
         for c = 1:size(group.avg.dPrime,1)
            % fit group-average
            polyparams = polyfit(data.ecc,group.avg.dPrime(c,:),data.model.polyOrder);
            fitVal     = polyval(polyparams,data.model.ecc);

            % draw fit
            plot(data.model.ecc,fitVal,'-','linewidth',4,'color',colors(c,:)); hold on
            
            % add area for peak eccentricity
            forCI = linspace(group.err.peakEcc(c,1),group.err.peakEcc(c,2),1e2);
            h = area(forCI,polyval(polyparams,forCI));
            set(h,'FaceColor',colors(c,:));
            set(h,'FaceAlpha',0.1);
            set(h,'EdgeAlpha',0);
            
            % raw data
            leg(c) = plot(data.ecc,group.avg.dPrime(c,:),'.','markersize',20,'color',colors(c,:)); hold on

            % errobars
            errorbar(data.ecc,group.avg.dPrime(c,:),group.err.dPrime(c,:),'color',colors(c,:),'linestyle','none','capsize',0,'linewidth',1.5);
         end
         % pretty up figure
         figureDefaults
         set(gca,'xtick',0:2:8,'xlim',[-0.3 8],'ytick',0:.5:2,'ylim',[0.4 2.2],'ticklength',[0.025 0.05]);
         ylabel('Performance (d^{\prime})','fontname','arial','fontsize',10);
      
         
      % reaction time
         subplot(2,1,2)
         
         for c = 1:size(group.avg.rt,1)
            % raw data
            leg(c) = plot(data.ecc,group.avg.rt(c,:),'.-','markersize',20,'color',colors(c,:),'linewidth',2); hold on

            % errobars
            errorbar(data.ecc,group.avg.rt(c,:),group.err.rt(c,:),'color',colors(c,:),'linestyle','none','capsize',0,'linewidth',1.5);
         end
         % pretty up figure
         figureDefaults
         set(gca,'xtick',0:2:8,'xlim',[-0.3 8],'ytick',0.1:0.1:1,'ylim',[0.1 0.22],'ticklength',[0.025 0.05]);
         xlabel('Eccentricity (visual angle)','fontname','arial','fontsize',10);
         ylabel('RT (s)','fontname','arial','fontsize',10);
         legend(leg,{'Neutral' 'Central'},'fontname','arial','fontsize',8,'location','southwest');
