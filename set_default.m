   
function defaults = set_default(function_name)
% create a structure with the default parameters of various 
% functions used by the boot_tide package

    switch function_name
        case 'boot'
        
            % set defaults  bootstrap parameters
            defaults.nboot       = 10^3;        % number of bootstrap resample to generate
            defaults.seed        = 2*floor(rand(2,1)*(10^8)) +1; % pseudo-numebr generator seed: generate an odd number as seed  
            defaults.method      = 'MBB';       % bootstrap type: MBB or SPB
            defaults.circular    = false;       % use circular MBB or not
            defaults.lblock_dist = 'geom';      % default distribution (MATLAB name) for simulating the clock length
            defaults.lblock_par  = [];          % block length distribution parameter 
            defaults.lblock      = 31*24;       % MBB block length [h]

        case {'out','output'}
            defaults.options   = NaN;
            defaults.yhat      = [];
            defaults.yres      = [];   

        otherwise
            error('set_default: unrecognized function_name');      
    end
end