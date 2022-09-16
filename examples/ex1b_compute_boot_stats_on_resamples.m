% ex1b_compute_boot_stats_on_resamples.m
%
% Use boot_tide.m for generating MBB water level resamples
% at two stations and compute bootstrap plug-in statistics with boot_tide_param.m.
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 
clc; clear

% load one year of hourly data at two stations:
load_two_stations

% MBB method
disp('>> apply boot_tide MBB without options (default values)')
mbb_out = boot_tide(t(:,1),h,'yhat',hhat); 

% water level resamples
wl_resamples = mbb_out.y_boot;

% compute the bootstrap plug-in estimates
tic; [wl_mean,wl_std,wl_ci] = boot_tide_param(wl_resamples); tc =toc;

% print the computing time
disp(['computing boot statistics from ' num2str(mbb_out.options.nboot) ' water level resamples took ' num2str(tc) ' sec.'])


    
% construct graphics for the mean water level simulated over a short period (one week)
% (just to check the series properties)

fig = figure('visible','off'); % def a figure obj

         t0 = datenum(2000,09,01,0,0,0);  % date start
         te = datenum(2000,09,07,23,0,0); % date end
        Tgr = ((t0*24:te*24)/24)';          
[tgr,it,iT] = intersect(roundn(t,-3),roundn(Tgr,-3));
       indt = 1:48:length(tgr);
       
       
% create graphics for a random week
for s = 1 : 2
    subplot(2,1,s); 
    hold on
    
     yorig = h(it,s);
     yboot = wl_mean(it,s);
     plot(tgr,yorig);
     plot(tgr,yboot);


         
   set(gca,'xtick',tgr(indt),'xticklabel',datestr(tgr(indt)))
   ylabel('water level [m]')
   grid on
   title(stations{s})

end       
legend('original','boot avg')      
print(fig,'-dpdf',['figs/gr1_mean_boot_wat_lev'],'-fillpage') 
close all


% construct graphics for the width of water level ci estimated over a short period (one week)
fig = figure('visible','off'); % def a figure obj


% create graphics for a random week
for s = 1 : 2
    subplot(2,1,s); 
    hold on
    
     yorig     = h(it,s);
     yboot_ciw = wl_ci(it,2,s) - wl_ci(it,1,s);
     plot(tgr,yboot_ciw);
         
   set(gca,'xtick',tgr(indt),'xticklabel',datestr(tgr(indt)))
   ylabel('water level CI width [m]')
   grid on
   title(stations{s})

end       
print(fig,'-dpdf',['figs/gr1_wat_lev_ci'],'-fillpage') 
close all


