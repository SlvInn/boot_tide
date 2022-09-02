   
function defaults = set_default(function_name)
% create a structure with the default parameters of various 
% functions used by the boot_tide package

    switch function_name

        case 'data'

            % define a data structure
            defaults.t      = [];   % time
            defaults.y      = [];   % observations
            defaults.yhat   = [];  % reconstructions
            defaults.yres   = [];  % residuals
            defaults.valid  = NaN; % index of valid data (non NaN time, observation, and residual values)
            defaults.tresol = 1;   % data temporal resolution [h]
    
        case 'boot'
        
            % set defaults  bootstrap parameters
            defaults.nboot        = 10^3;        % number of bootstrap resample to generate
            defaults.seed         = 2*floor(rand(2,1)*(10^8)) +1; % pseudo-numebr generator seed: generate an odd number as seed  
            defaults.method       = 'mbb';       % bootstrap type: MBB or SPB
            defaults.circular     = false;       % use circular MBB or not
            defaults.noise        = 'fft';       % method to estimate the residual spectrum: fft, fftnoise, skfft, sknoise, cpfft, cpnoise, dllftt, or dllnoise
            
            defaults.lblock_dist  = 'geom';      % default distribution (MATLAB name) for simulating the clock length
            defaults.lblock_par   = [];          % block length distribution parameter 
            defaults.lblock       = 31*24;       % MBB block length [h]
            defaults.block_output = false;       % if true return the MBB matrices used for construction the resamples (mbb block sampling)
            
            defaults.tide_model = [];  % function handle for eatimating the tidal parameters of the y resamples (y_boot)
            defaults.theta      = [];  % vector of parameter estimated on the observed sample [i.e., theta = tide_model(y)]

        case {'out','output'}
            
            defaults.yhat      = [];   % reconstructed / predicted water levels at times t
            defaults.yres      = [];   % residuals between obs and predicted water levels: yres = y - yhat
            defaults.options   = NaN;  % options used by the boot_tide.m call


        case 'boot_param'

            defaults.circular   = [];          % vector indicating which compontent of theta are circural variables (e.g., phases)  
            defaults.circ_units = 'rad';       % units of circular parameters: 'rad', 'radians', 'deg', or 'degrees'
            defaults.ci_method = 'percentile' % string indicating the confidence interval method to use: 'percentile', 'gaussian', 'accelerated', 'bca'
            defaults.alpha      = 0.05;        % confidence level for constructing CI
 
        otherwise
            error('set_default: unrecognized function_name');      
    end
end