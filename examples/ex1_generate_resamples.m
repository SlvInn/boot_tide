% ex1_generate_resample.m
% matlab -nodisplay -batch ex1_generate_resample.m
% 
% Test the boot_tide.m fonction for generating MBB and SPB residual resamples
% for on exemple of water level reconstructions.
% Note: aux() is a function container used to store miscellaneaous (sub)-functions in aux.m 
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 
clc
clear
addpath('../')

% fix the period
t0 = datenum(2000,01,01,0,0,0);  % date start
te = datenum(2000,12,31,23,0,0); % date end

% load data at 2 stations 
stations = {'Halifax', 'St.Francois'};
nS = length(stations);

% load the observed water levels
[h(:,1),t] = aux('read_hourly_data','./data/hobs_HAL_1999_2009.mat',t0,te); % water level and time @ Halifax
[h(:,2),~] = aux('read_hourly_data','./data/hobs_STF_1999_2009.mat',t0,te); % water level @ St.Francois

% load some reconstructions at these two stations (from HA)
[hhat(:,1),t] = aux('read_hourly_data','./data/hhat_HAL_2000.mat',t0,te); % estimated water level and time @ Halifax
[hhat(:,2),~] = aux('read_hourly_data','./data/hhat_STF_2000.mat',t0,te); % estimated water level @ St.Francois


# use a very low nboot value to have faster calculations
nboot = 1000;

% apply boot_tide MBB with defaults values, except nboot 
% mbb_dst = {'unif' 'pois' 'geom' 'fix'};
% for dd = 1 : 4
%     disp(['MBB ' mbb_dst{dd}])
%     out_mbb.(mbb_dst{dd}) = boot_tide(t,h,'yhat',hhat,'nboot',nboot,'lblock_dist',mbb_dst{dd});
% end

% % apply boot_tide SPB with defaults values, except nboot 
% spb_mtd = {'fft','skfft', 'cpfft', 'dllftt'};
% for dd = 1 : 4
%     disp(['SPB ' spb_mtd{dd}])
%     out_spb.(spb_mtd{dd}) = boot_tide(t,h,'yhat',hhat,'nboot',nboot,'method','spb','noise',spb_mtd{dd});
% end


% save the test output
% var2save = {'out_mbb','out_spb'};
% save('./data/test_default_boot_tide.mat',var2save{:})
load('./data/test_default_boot_tide.mat')


% SILVIA: sono qui >> devo testare from_mbb_to_resamples.m
% POI:
% - fare una funzione boot per il calcolo dei parametri (plug-in) > lin vs circ
% - fare una funzione che sa applicare qualsiasi metodo di stima ? 
% - fare un template per l'uso del package
% - compara con NS Tide
% - fai grafici per le ricostruzioni, i bloccchi, compara i metodi SPB, etc

