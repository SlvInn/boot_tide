   
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
        
            % set default bootstrap parameters
            defaults.nboot        = 10^4;        % number of bootstrap resamples to generate
            defaults.seed         = 2*floor(rand(2,1)*(10^8)) +1; % pseudo-number generator seed: generate an odd number as seed  
            defaults.method       = 'mbb';       % bootstrap type: MBB or SPB
            defaults.circular     = true;        % use circular MBB or not
            defaults.noise        = 'fft';       % method to estimate the residual spectrum: fft, fftnoise, skfft, sknoise, cpfft, cpnoise, dllftt, or dllnoise
            
            defaults.lblock_dist  = 'geom';      % default distribution (MATLAB name) for simulating the block length
            defaults.lblock_par   = [];          % block length distribution parameter 
            defaults.lblock       = 31*24;       % MBB block length [h]
            defaults.block_output = false;       % if true return the MBB matrices used for constructing the resamples (mbb block sampling)
            
            defaults.tide_model = [];  % function handle for estimating the tidal parameters on each y resample 
            defaults.theta      = [];  % vector of parameter estimated on the observed sample [i.e., theta = tide_model(y)]
            defaults.iboot      = [];  % lt x nboot matrix with temporal index to use to construct the resamples
            

        case {'out','output'}
            
            defaults.yhat      = [];   % reconstructed / predicted water levels at times t
            defaults.yres      = [];   % residuals between obs and predicted water levels: yres = y - yhat
            defaults.options   = NaN;  % options used by the boot_tide.m call


        case 'boot_param'

            defaults.circular   = [];           % vector indicating which component of theta are circural variables (e.g., phases)  
            defaults.circ_units = 'rad';        % units of circular parameters: 'rad', 'radians', 'deg', or 'degrees'
            defaults.ci         = 'percentile'; % string indicating the confidence interval method to use: 'percentile', 'gaussian', 'studentized', 'accelerated', 'bca'
            defaults.alpha      = 0.05;         % confidence level for constructing CI
            defaults.theta0      = [];           % vector of paramaters computed on the ortiginal samples (used for studentizzed boot CI)

        otherwise
            error('set_default: unrecognized function_name');      
    end
end