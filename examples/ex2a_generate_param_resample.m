% ex2_generate_param_resample.m
% 
% Test the boot_tide.m fonction for generating resamples of the model 
% parameters using the MBB. Then compute the parameter bootstrap (plug-in) 
% estimates and CI.
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 
clc; clear

% load one year of hourly data at one station:
load_one_station



%% define a model for computing the bootstrap parameter replicates
% This allows to estimate the models on all the resamples without 
% calling the original tide analysis package (much faster since we avoid most of the 
% time consuming steps)


% use some quantities estimated on the original water level sample 
w = sqrt(cf.W);     % IRLS weight
theta = cf.M;       % HA parameters (sin-cos amplitudes and the intercept)
B1 = [B ones(size(B,1),1)]; % add a column for the intercept
  
% function handle to be passed to tide_model 
tidemdl = @(y) (w(~isnan(y)).*B1(~isnan(y),:)) \ (w(~isnan(y)).*y(~isnan(y))); 
 

%% apply the MBB

% define the boot parameters
nboot = 10^2;
boot_opt = {'yhat',hhat,'nboot',nboot,'lblock_dist','geom','tide_model',tidemdl,'theta',theta};
% use a low value of nboot to reduce the computing time

% construct the resamples
tic; mbb_out = boot_tide(t,h,boot_opt{:}); 
tc_mbb_wtm = toc; % stock the computing time 

% constructed resamples in mbb_out.theta_boot (HA regression coeffs):
disp('constructed resamples of parameters in mbb_out.theta_boot (HA regression coeffs):')
disp(mbb_out)

% print the computing time
disp(['generating ' num2str(nboot) ' HA regress coeff resamples took ' num2str(tc_mbb_wtm) ' sec.'])



% estimate the boot estimators (stats and CI) for the raw parameters (HA regression coeffs)
tic; 
beta_boot = mbb_out.theta_boot;
[m_beta,s_beta] = boot_tide_param(beta_boot); 


% estimate the tidal amp and phases from each beta replicate 
for b = 1 : nboot
    [mbb_out.a_boot(:,b),mbb_out.g_boot(:,b)] = from_bsincos_to_amppha(beta_boot(:,b));
end

% estimate the boot estimators  for the amplitudes
[m_amp,s_amp, ci_amp] = boot_tide_param(mbb_out.a_boot);


% create a boolean vector to set 'circular' to true for each phase param
circ = mbb_out.g_boot(:,1).*0+1;

% estimate the boot stats for the phases
[m_pha,s_pha,ci_pha] = boot_tide_param(mbb_out.g_boot,'circular',circ,'circ_units','deg');
tc_boot_ci = toc; % stock the computing time 

% print the computing time
disp(['estimating amplitude and phase CIs from the ' num2str(nboot) ' HA regress coeff resamples took ' num2str(tc_boot_ci) ' sec.'])


% check the computated bootstrap statistics (estimators) and CI on graphics
% ex2b_generate_param_resample_fig.m