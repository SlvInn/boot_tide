% ex3a_compare_boot_ci.m
%
% Compute the tidal amplitude and phase boot estimators and CI
% with the percentile, gaussian, and studentized boot methods 
% and compare them graphically.
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 
clc; clear

% load one year of hourly data at two stations:
load_two_stations

% MBB method
disp('>> apply boot_tide MBB without options (default values)')
mbb_out = boot_tide(t(:,1),h,'yhat',hhat); 


% use some quantities estimated on the original water level sample 
theta = [coef.hal.M, coef.stf.M];     % HA parameters (sin-cos amplitudes and the intercept)
B1 = [coef.hal.B ones(size(coef.hal.B,1),1)]; % use the same regressor for both stations

% function handle to be passed to tide_model 
tidemdl = @(y) B1(~isnan(y),:) \ y(~isnan(y)); 
% NOTE: consider a simple OLS regression since the same tidal_model is applied to both stations


% define the boot parameters
nboot    = 40; %10^2;
boot_opt = {'yhat',hhat,'nboot',nboot,'lblock_dist','geom','tide_model',tidemdl,'theta',theta};
mbb_out  = boot_tide(t(:,1),h,boot_opt{:}); % construct the resamples


% at each station, estimate the tidal amp and phases from each beta replicate 
for s = 1 : ns
    for b = 1 : nboot
        [mbb_out.a_boot(:,b,s),mbb_out.g_boot(:,b,s)] = from_bsincos_to_amppha(mbb_out.theta_boot(:,b,s));
    end
end


% estimate the boot estimators  for the amplitudes with different CI formulas
ci_methods = {'percentile','gaussian','studentized' };
A0 = [coef.hal.A,coef.stf.A]; % theta0 needed for the studentized boot

for s = 1 : ns
    for c = 1 : 3
        mtd = ci_methods{c};
        [m_amp.(mtd),s_amp.(mtd), ci_amp.(mtd)] = boot_tide_param(mbb_out.a_boot,'ci',mtd,'theta0',A0);
    
        % NOTE: the gaussian and studentized method allow for ci_amp<0 
        % since the both use a gaussian distribution for the amplitudes
        % >> put to NaN negative ci_amp values:
        ci_amp.(mtd)(ci_amp.(mtd)<0)= nan;
    
    end
end



% construct a graphic for the log of the amplitude CI bounds (in cm)

% boot amplitudes
ticks = false;
col   = [ 27,158,119; 217,95,2;0 0 0; ]./255; %  150 150 150;
bound_names = {'CI lower bound','CI upper bound'};
for s = 1: ns
    fig = figure('visible','off'); % def a figure obj

    for c = 1 : 3
        mtd = ci_methods{c};


        % amplitudes boot CIs
        ciA = log(ci_amp.(mtd)*100);

        

        for bound = 1 : 2

            if c==3 && bound==2
                ticks = true;
            end

            subplot(2,1,2-bound+1); hold on
            plot_amplitudes(ciA(:,bound,s),coef.(stn{s}).name,coef.(stn{s}).aux.frq,...
                'ticks',ticks,'staircase',true,'-','color',col(c,:),'linewidth',1);
            
            % set the axes formats
            if c==3 
                ytk = [0.0005 0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.2 0.5 1]*100;
                set(gca,'ytick',log(ytk),'yticklabel',ytk)
                ylim(log([0.005 220]))

                grid on 
                ylabel('A_k [cm]')
                legend(ci_methods,'location','NE')

            end
            title([bound_names{bound}])
        end
    end

    print(fig,'-dpdf',['figs/gr3_boot_CIs_' stations{s}],'-fillpage') 
close all

end


