% ex1_generate_wl_resample.m
% matlab -nodisplay -batch ex1_generate_resample
% 
% Test the boot_tide.m fonction for generating MBB and SPB residual resamples
% for on exemple of water level reconstructions.
% Note: aux() is a function container used to store miscellaneaous (sub)-functions in aux.m 
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 
clc; clear

% load the data using the r1_load_stations.m script:
r1_load_stations

% use a very low nboot value to have faster calculations
nboot = 10^4; %10^2; %

% % % apply boot_tide MBB with defaults values, except nboot 
mbb_dst = {'unif' 'pois' 'geom' 'fix'};
for dd = 1 : 4
    disp(['MBB ' mbb_dst{dd}])
    tic
    out_mbb.(mbb_dst{dd}) = boot_tide(t,h,'yhat',hhat,'nboot',nboot,'lblock_dist',mbb_dst{dd});
    comp_time.(mbb_dst{dd}) = toc; % computing time
end

% % apply boot_tide SPB with defaults values, except nboot 
spb_mtd = {'fft','skfft', 'cpfft', 'dllftt'};
for dd = 1 : 4
    disp(['SPB ' spb_mtd{dd}])
    tic
    out_spb.(spb_mtd{dd}) = boot_tide(t,h,'yhat',hhat,'nboot',nboot,'method','spb','noise',spb_mtd{dd});
    comp_time.(spb_mtd{dd}) = toc;
end


% % save the test output
var2save = {'out_mbb','out_spb','comp_time'};
save('./data/wl_default_resamples.mat',var2save{:})


% load the estimated data
% load('./data/wl_default_resamples.mat')




% mbb_dst = {'unif' 'pois' 'geom' 'fix'};
% spb_mtd = {'fft','skfft', 'cpfft', 'dllftt'};
% for dd = 1 : 4
%     out_mbb.(mbb_dst{dd}).y_boot = out_mbb.(mbb_dst{dd}).theta_boot; 
%     out_spb.(spb_mtd{dd}).y_boot = out_spb.(spb_mtd{dd}).theta_boot; 
% end
% save('./data/wl_default_resamples_100b.mat',var2save{:})
% load('./data/wl_default_resamples_100b.mat')






