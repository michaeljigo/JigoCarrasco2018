% Purpose:  Recreate Figure 4. 
%           Cueing effects for exogenous and endogenous attention conditions.
%
% By:       Michael Jigo
% Edited:   07.01.21
%
% Output:   stats    structure containing results of statistical analyses

function stats = figure4
addpath(genpath('./helperfun'));

   %% Get subject-level results
      [data params]  = pack_experiment_results;
      attention_type = {'exo' 'endo'};
      for a = 1:numel(attention_type)
         idx            = find(ismember(params.attntype,attention_type{a})); 
         effect(:,a,:)  = data(idx).cueing_effect;
      end
      data = data(1);
      data.cueing_effect = effect;


   %% Compute group-average and error
      % cueing effect
         group.avg.effect  = squeeze(nanmean(data.cueing_effect,1));
         group.err.effect  = withinSubjError(data.cueing_effect,0);


   %% Perform stats
      % d-prime
         stats.effect      = simple_mixed_anova(data.cueing_effect,[],{'attention' 'ecc'});


   %% Plot
      colors = [0 0 0; 202 0 32]./255;
      figure('name','Cueing effects','position',[617 -8 336 553]);
      
      for c = 1:size(group.avg.effect,1)
         % line at 0
         line([-0.3 8],[0 0],'color',[0.5 0.5 0.5],'linestyle','--','linewidth',1.5);

         % fit group-average
         polyparams = polyfit(data.ecc,group.avg.effect(c,:),1);
         fitVal     = polyval(polyparams,data.model.ecc);

         % draw fit
         plot(data.model.ecc,fitVal,'-','linewidth',4,'color',colors(c,:)); hold on
            
         % raw data
         leg(c) = plot(data.ecc,group.avg.effect(c,:),'.','markersize',20,'color',colors(c,:)); hold on

         % errobars
         errorbar(data.ecc,group.avg.effect(c,:),group.err.effect(c,:),'color',colors(c,:),'linestyle','none','capsize',0,'linewidth',1.5);
      end
      % pretty up figure
      figureDefaults
      set(gca,'xtick',0:2:8,'xlim',[-0.3 8],'ytick',-0.3:0.3:0.9,'ylim',[-0.35 0.9],'ticklength',[0.025 0.05]);
      ylabel('Valid-Neutral (\Delta d^{\prime})','fontname','arial','fontsize',10);
      xlabel('Eccentricity (visual angle)','fontname','arial','fontsize',10);
      legend(leg,{'Central' 'Peripheral'},'fontname','arial','fontsize',8,'location','southeast');
