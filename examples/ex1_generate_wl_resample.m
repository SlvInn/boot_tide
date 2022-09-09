% ex1_generate_wl_resample.m
% 
% Test the boot_tide.m fonction for generating MBB and SPB residual resamples
% for on exemple of water level reconstructions. 
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 
clc; clear

% load the data using the r1_load_stations.m script:
r1_load_stations

% use a low nboot value to have faster calculations in tests
nboot = 1000; 



% compute bootstrap resamples
%

% % % apply boot_tide MBB with defaults values, except nboot and the
% % % distribution of the block length 
mbb_mtd = {'unif' 'pois' 'geom' 'fix'};
for dd = 1 : 4
    disp(['MBB ' mbb_mtd{dd}])
    tic
    out_mbb.(mbb_mtd{dd}) = boot_tide(t,h,'yhat',hhat,'nboot',nboot,'lblock_dist',mbb_mtd{dd});
    comp_time.(mbb_mtd{dd}) = toc; % computing time
end


% % % apply boot_tide SPB with defaults values, except nboot and the
% % % function used to simulate the residual resamples via their sprectrum
spb_mtd = {'fftnoise','sknoise', 'cpnoise', 'dllnoise'}; % ftt method for generating the SPB noise 
for dd = 1 : 4
    disp(['SPB ' spb_mtd{dd}])
    tic
    out_spb.(spb_mtd{dd}) = boot_tide(t,h,'yhat',hhat,'nboot',nboot,'method','spb','noise',spb_mtd{dd});
    comp_time.(spb_mtd{dd}) = toc;
end



% % save the test output
var2save = {'out_mbb','out_spb','comp_time'};
save('./data/wl_default_resamples.mat',var2save{:})
%}

% load the estimated HA and teh bootstrap resamples
load('./data/wl_default_resamples.mat')


% combine MBB and SPB structures
  out = out_mbb;
    f = fieldnames(out_spb);
for i = 1:length(f)
    out.(f{i}) = out_spb.(f{i});
end
mtd = fieldnames(out);
clear out_mbb out_spb


% compute some bootstrap statistics on the resamples
mtd  = fieldnames(out);
nMtd = length(mtd);

for i = 1:nMtd
    
    % compute the bootstrap plug-in estimates
    [mh,sh,h_ci] = boot_tide_param(out.(mtd{i}).y_boot);
    
    
     y_boot(:,i,:) = mh;          % avg of the water lev series reconstructed via the bootstraps
     ci_low(:,i,:) = h_ci(:,1,:); % lower CI limit for water level reconstructions
     ci_up(:,i,:)  = h_ci(:,2,:); % upper CI limit for water level reconstructions
     
end



% construct graphics for the mean water level simulated over a short period (one week)
fig = figure; % def a figure obj

         t0 = datenum(2000,09,01,0,0,0);  % date start
         te = datenum(2000,09,07,23,0,0); % date end
        Tgr = ((t0*24:te*24)/24)';          
[tgr,it,iT] = intersect(roundn(t,-3),roundn(Tgr,-3));
       indt = 1:48:length(tgr);

for i = 1:nMtd
         
    % create graphics for a random week
    for s = 1 : 2
        subplot(2,1,s); hold on
        ygr = y_boot(it,i,s);
        pl(s,i)= plot(tgr,ygr);
        
        
        if i == nMtd         
           set(gca,'xtick',tgr(indt),'xticklabel',datestr(tgr(indt)))
           grid on
           title(stations{s})
        end
    end
end
legend(pl(1,:),mtd{:}); % add the method names in the legend
print(fig,'-dpdf',[FldOutGr '\gr1_wat_lev_mean'],'-fillpage') 


% construct graphics for the width of water level ci estimated over a short period (one week)
fig = figure('visible','off'); % def a figure obj

for i = 1:nMtd
         
     % create graphics for a random week
    for s = 1 : 2
        subplot(2,1,s); hold on
        ygr = ci_up(it,i,s) - ci_low(it,i,s);
        pl(s,i)= plot(tgr,ygr);
        
        
        if i == nMtd         
           set(gca,'xtick',tgr(indt),'xticklabel',datestr(tgr(indt)))
           grid on
           title(stations{s})
        end
    end
end
legend(pl(1,:),mtd{:}); % add the method names in the legend
print(fig,'-dpdf',[FldOutGr '\gr1_wat_lev_ci'],'-fillpage') 
close all