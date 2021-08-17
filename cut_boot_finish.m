
function coef_boot = cut_boot_finish(coef_boot,opt)
% re-order the fields of the coed_boot structure

if opt.twodim
    fldord = {'name'; 'M'; 'M_rlzn';  'Lsmaj'; 'Lsmaj_rlzn'; ...  
              'Lsmin'; 'Lsmin_rlzn'; 'theta';  'theta_rlzn';  'g';  'g_rlzn'; ...
              'umean'; 'umean_rlzn'; 'vmean'; 'vmean_rlzn'; 'u';'v';'W'; 'B'; 'Std'; 'boot'; 'boot_opt'};
else
    fldord = {'name';  'M'; 'M_rlzn'; 'A'; 'A_rlzn'; ...  
               'g';  'g_rlzn';  'mean'; 'mean_rlzn'; ...
               'u';'W'; 'B'; 'Std'; 'boot';'boot_opt'}; 
end
if ~opt.notrend
    if opt.twodim
        fldord{end+1} = 'uslope';
        fldord{end+1} = 'uslope_rlzn';
        fldord{end+1} = 'vslope';
        fldord{end+1} = 'vslope_rlzn';
    else
        fldord{end+1} = 'slope';
        fldord{end+1} = 'slope_rlzn';
    end
end
   
coef_boot = orderfields(coef_boot,fldord);
      
end
