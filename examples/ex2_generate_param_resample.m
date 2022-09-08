% ex2_generate_param_resample.m
% 
% Test the boot_tide.m fonction for generating resamples of the model 
% parameters using the MBB and SPB. Then compute the parameter bootstrap (plug-in) 
% estimates and CI.
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 
clc; clear

% load the data using the r1_load_stations.m script:
r2_load_halifax_station

% use a low nboot value to have faster calculations in tests
nboot = 10^3; 


% some quantities estimated on the original water level sample 
      w = sqrt(cf.W); % IRLS weight
  theta = cf.M;       % HA parameters (sin-cos amplitudes and the intercept)
     B1 = [B ones(size(B,1),1)]; % add a column for the intercept
     
% function/model to be used to compute the bootstrap parameter replicates    
tidemdl = @(y) (w(~isnan(y)).*B1(~isnan(y),:)) \ (w(~isnan(y)).*y(~isnan(y))); 

% construct the resamples
%{
% put common parameters together
boot_opt = {'yhat',hhat,'nboot',nboot,'tide_model',tidemdl,'theta',theta};

% % % apply boot_tide MBB with defaults values, except nboot and the
% % % distribution of the block length 
mbb_mtd = {'unif' 'pois' 'geom' 'fix'};
for dd = 1 : 4 % for each possible distribution simulating the block length
    disp(['MBB ' mbb_mtd{dd}])
    tic
    out.(['MBB_' mbb_mtd{dd}]) = boot_tide(t,h,'lblock_dist',mbb_mtd{dd},boot_opt{:});
    comp_time.(mbb_mtd{dd}) = toc; % computing time
end

% % % apply boot_tide SPB with defaults values, except nboot and the
% % % function used to simulate the residual resamples via their sprectrum
spb_mtd = {'fftnoise','sknoise', 'cpnoise', 'dllnoise'}; % ftt method for generating the SPB noise
for dd = 1 : 4
    disp(['SPB ' spb_mtd{dd}])
    tic
    out.(['SPB_' spb_mtd{dd}]) = boot_tide(t,h,'method','spb','noise',spb_mtd{dd},boot_opt{:});
    comp_time.(spb_mtd{dd}) = toc; % computing time
end



% % save the output
var2save = {'out','comp_time'};
save('./data/ha_param_resamples.mat',var2save{:})
%}

% load the estimated HA and the bootstrap parameter replicates
load('./data/ha_param_resamples.mat')


mtd  = fieldnames(out);
nMtd = length(mtd);

for i = 1:nMtd
    
    % estimate the boot parameters and CI for the raw parameters
    th_b = out.(mtd{i}).theta_boot;
    [mM, sM]   = boot_tide_param(th_b); 
    theta_boot.(mtd{i}).avg.M = mM;
    theta_boot.(mtd{i}).std.M = sM;
    
    % estimate the amp and phases from each param replicate 
    for b = 1 : nboot
      [out.(mtd{i}).a_boot(:,b),out.(mtd{i}).g_boot(:,b)] = aux('from_bsincos_to_amppha',th_b(:,b));
    end
    
    % estimate the boot stats for the amplitudes
    [mA,sA, ciA] = boot_tide_param(out.(mtd{i}).a_boot);
    theta_boot.(mtd{i}).avg.A     = mA;
    theta_boot.(mtd{i}).std.A = sA;
    theta_boot.(mtd{i}).ci.A  = ciA;
    
    
    % create a boolean vector to set circular to true for each phase param
    circ = out.(mtd{i}).g_boot(:,1).*0+1;
    
    % estimate the boot stats for the phases
    [mg,sg,cig] = boot_tide_param(out.(mtd{i}).g_boot,'circular',circ,'circ_units','deg');
    theta_boot.(mtd{i}).avg.g = mg;
    theta_boot.(mtd{i}).std.g = sg;
    theta_boot.(mtd{i}).ci.g  = cig;
end

% check the computated bootstrap statistics (estimators) and CI on graphics
fig   = figure;
ticks = false;
for i = 1:nMtd

  
    % construct a graphic for the log of the amplitudes (in cm)
    subplot(2,1,1); hold on
    mA = log(theta_boot.(mtd{i}).avg.A*100);
    pl(i) = aux('plot_amplitudes',mA,cf.name,cf.aux.frq,ticks);
    
    % add a ref and set the axes formats
    if i==nMtd 
       mA = log(cf.A*100);
          aux('plot_amplitudes',mA,cf.name,cf.aux.frq,true,':k'); % UTide Ak
       ylabel('log(A_k)')
        title('Tidal amplitude')
    end
 
     % construct a graphic for the log of the amplitude stds (in cm)    
    subplot(2,1,2); hold on
    mA = log(theta_boot.(mtd{i}).std.A*100);
    aux('plot_amplitudes',mA,cf.name,cf.aux.frq,ticks);
    
    % add a ref and set the axes formats
    if i==nMtd
       mA = log(cf.Std.A*100);
          aux('plot_amplitudes',mA,cf.name,cf.aux.frq,true,':k');% UTide std(Ak)
       ylabel('log(\sigma_{A_k})')
        title('Tidal amplitude std')
    end
end
legend(pl(1,:),mtd{:}); % add the method names in the legend

