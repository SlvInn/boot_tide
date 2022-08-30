function coef_boot = cut_boot_stats(coef_boot,n_boot,opt,nNR,nR)
% compute the component-wise mean of M and the corresponding 
% amplitude, phase, mean, and trend parameters; then compute the std 
% of linear parameters (amplitudes,means and trends) and the circular
% variance of phases


    % compute the component-wise mean of M 
    M = coef_boot.M_rlzn;
    nMr =  size(M,1);
    coef_boot.M = nan(nMr,1);
    for r = 1 : nMr
        Mr = M(r,:);
   coef_boot.M(r) = mean(Mr(~isnan(Mr)));
    end

    % compute the amplitude, phase, mean, and trend parameters corresponding to the mean M
       mcoef_b = coef_boot;
        F_rlzn = fields(mcoef_b);
         for f = 1 : length(F_rlzn)
            if contains(F_rlzn{f},'_rlzn')
                 mcoef_b.(F_rlzn{f}) = [];
            end
        end
        mcoef_b.M_rlzn = coef_boot.M;
        mcoef_b = cut_from_mboot_to_param(mcoef_b,1,opt,nNR,nR);


% compute the parameters corresponding to the mean M
coef_boot.g = mcoef_b.g_rlzn;
if ~opt.twodim                      
    coef_boot.A    = mcoef_b.A_rlzn;
    coef_boot.mean = mcoef_b.mean_rlzn;
    if ~opt.notrend
        coef_boot.slope = mcoef_b.slope_rlzn;
    end        
else
    coef_boot.Lsmaj = mcoef_b.Lsmaj_rlzn;
    coef_boot.Lsmin = mcoef_b.Lsmin_rlzn;
    coef_boot.theta = mcoef_b.theta_rlzn;
    coef_boot.umean = mcoef_b.umean_rlzn;
    coef_boot.vmean = mcoef_b.vmean_rlzn;
    if ~opt.notrend
        coef_boot.uslope = mcoef_b.uslope_rlzn;
        coef_boot.vslope = mcoef_b.vslope_rlzn;
    end
end

% compute the variability of the estimates
nC  = numel(coef_boot.name);
mad = @(x) median(abs(x-median(x,2)),2)./0.6745; 

coef_boot.Std.M  = mad(coef_boot.M_rlzn);

  d2r = (pi./180);
for c = 1 : nC
    coef_boot.Std.g(c,1) = circ_var(d2r*coef_boot.g_rlzn(c,:)'); 
end

if ~opt.twodim                      
    coef_boot.Std.A    = mad(coef_boot.A_rlzn);
    coef_boot.Std.mean = mad(coef_boot.mean_rlzn);
    if ~opt.notrend
        coef_boot.Std.slope = mad(coef_boot.slope_rlzn);
    end        
else
    for c = 1 : nC
        coef_boot.Std.theta(c,1) = circ_var(d2r*coef_boot.theta_rlzn(c,:)');
    end
    coef_boot.Std.Lsmaj = mad(coef_boot.Lsmaj_rlzn);
    coef_boot.Std.Lsmin = mad(coef_boot.Lsmin_rlzn);
    
    coef_boot.Std.umean = mad(coef_boot.umean_rlzn);
    coef_boot.Std.vmean = mad(coef_boot.vmean_rlzn);
    if ~opt.notrend
        coef_boot.Std.uslope = mad(coef_boot.uslope_rlzn);
        coef_boot.Std.vslope = mad(coef_boot.vslope_rlzn);
    end
end



end