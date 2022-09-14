% ex2b_generate_param_resample_fig.m 
%
% Compute the tidal amplitude and phase boot estimators and CI
% and plot them on a graphic.
%
% silvia.innocenti@ec.gc.ca, 2022/09/14. 
clc; clear

% load one year of hourly data at one station:
load_one_station



% use some quantities estimated on the original water level sample 
w = sqrt(cf.W);     % IRLS weight
theta = cf.M;       % HA parameters (sin-cos amplitudes and the intercept)
B1 = [B ones(size(B,1),1)]; % add a column for the intercept
  
% function handle to be passed to tide_model 
tidemdl = @(y) (w(~isnan(y)).*B1(~isnan(y),:)) \ (w(~isnan(y)).*y(~isnan(y))); 
 

% define the boot parameters
nboot    = 10^2;
boot_opt = {'yhat',hhat,'nboot',nboot,'lblock_dist','geom','tide_model',tidemdl,'theta',theta};
mbb_out  = boot_tide(t,h,boot_opt{:}); % construct the resamples

% estimate the tidal amp and phases from each beta replicate 
for b = 1 : nboot
    [mbb_out.a_boot(:,b),mbb_out.g_boot(:,b)] = from_bsincos_to_amppha(mbb_out.theta_boot(:,b));
end

% create a boolean vector to set 'circular' to true for each phase param
ag    = [mbb_out.a_boot; mbb_out.g_boot]; % amplitudes and phases
ncomp = length(mbb_out.a_boot(:,1));     % No. of tidal constituents
circ  = [mbb_out.a_boot(:,1).*0; mbb_out.g_boot(:,1).*0+1]; % boolean vector for 'circular' and 'linear' parameters

% estimate the boot estimators for the amplitudes and phases 
[m_ag,s_ag,ci_ag] = boot_tide_param(ag,'circular',circ,'circ_units','deg');



% check the computated bootstrap statistics and CI 
fig = figure('visible','off'); % def a figure obj


% construct a graphic for the log of the amplitudes (in cm)
subplot(2,1,1); hold on

% boot amplitudes
mA = log(m_ag(1:ncomp)*100);
pl(1) = plot_amplitudes(mA,cf.name,cf.aux.frq,'ticks',false,'staircase',true,'-k','linewidth',2);

% amplitudes boot CIs
ciA = log(ci_ag(1:ncomp,:)*100);
pl(2) = plot_amplitudes(ciA(:,1),cf.name,cf.aux.frq,'ticks',false,'staircase',true,'-k','linewidth',.5);
plot_amplitudes(ciA(:,2),cf.name,cf.aux.frq,'ticks',false,'staircase',true,'-k','linewidth',.5);

% add a ref 
mA = log(cf.A*100);
pl(3) = plot_amplitudes(mA,cf.name,cf.aux.frq,'ticks',true,'staircase',false,'xr'); % UTide estimation


% set the axes formats
ytk = [0.0005 0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.2 0.5 1]*100;
set(gca,'ytick',log(ytk),'yticklabel',ytk)
ylim(log([0.001 120]))

% text 
legend(pl(:),'mbb','mbb CI','UTide','location','S')
ylabel('A_k [cm]')
title('Tidal amplitude')

% construct a graphic for the phases (in degrees)
subplot(2,1,2); hold on

% boot phases
mg = m_ag(ncomp+1:end);
pl(1) = plot_amplitudes(mg,cf.name,cf.aux.frq,'ticks',false,'staircase',true,'-k','linewidth',2);

% amplitudes boot CIs
cig = ci_ag(ncomp+1:end,:);
pl(2) = plot_amplitudes(cig(:,1),cf.name,cf.aux.frq,'ticks',false,'staircase',true,'-k','linewidth',.5);
pl(3) = plot_amplitudes(cig(:,2),cf.name,cf.aux.frq,'ticks',false,'staircase',true,'-k','linewidth',.5);

% add a ref 
pl(4) = plot_amplitudes(cf.g,cf.name,cf.aux.frq,'ticks',false,'staircase',false,'xr'); % UTide estimation


% set the axes formats
ytk = (1:45:360);
set(gca,'ytick',ytk,'yticklabel',ytk)


% text 
% legend('','location','N')
ylabel('g_k [cm]')
title('Tidal phases')



print(fig,'-dpdf',['figs/gr2_boot_plugin_estimators'],'-fillpage') 
close all

