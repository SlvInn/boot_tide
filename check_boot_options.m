function opt = check_boot_options(data,opt)
% Apply basic checks on the entered bootstrap options
%
% INPUT:
% data - structure of processed data series   >> See help check_data
% opt  - structure with the bootstrap options >> See help boot_tide.m for a description of the opt structure fields 
% 
% OUTPUT:
% opt  - structure with the bootstrap options >> See help boot_tide.m
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2020/04 - 2022/08
 
    % number of bootstrap replicates/resamples
    assert(floor(opt.nboot)==opt.nboot,'nboot must be an integer')

    % random seed for resamplers
    assert(numel(opt.seed)==2,'the bootstrap random seed must be a 2 element vector')
    
    % bootstrap method: Moving Block Bootstrap (MBB) or Semi-Parametric Bootstrap (SPB)
    opt.method = lower(opt.method);
    assert(any(strcmp(opt.method,{'mbb','spb',''})),'not valid bootstrap method')

    
    if strcmp(opt.method,'mbb') % check options for the MBB
        opt = check_mbb_opt(data,opt);

        if opt.block_output>1
            warning('boot_tide:check_boot_options >> block_output>1 set to false')
            opt.block_output = logical(opt.block_output==true || opt.block_output<=1);
        end
        

    elseif strcmp(opt.method,'spb') % check options for SPB
        opt = check_spb_opt(opt); 
        
    else % check the iboot option and related ones
        opt = check_knownit_opt(opt); 

    end    

    % check that tide_model is a function handle or empty
    if ~isempty((opt.tide_model))
        assert(strcmpi(class(opt.tide_model),'function_handle'),'tide_model must be a function handle')
        assert(nargin(opt.tide_model)==1,'tide_model can only take 1 input argument (observation vector)')
        assert(~isempty((opt.theta)), 'theta must be provided for non empty tide_model')
        % opt.theta = opt.theta(:); %don't make theta being a col vector since we may have 2 stations
        
        ns = size(data.y,2);
        if ns>1
            warning(['applying the same tide_model to ' num2str(ns) ' spatial locations'])
        end
    else
        opt = rmfield(opt,{'tide_model','theta'});     
    end

end


function opt = check_mbb_opt(data,opt)

    % MBB circular option (avoid truncation at the series end)
    assert(numel(opt.circular)==1,'circular must be scalar (boolean or 0-1)')
    assert(islogical(opt.circular) || opt.circular==0 || opt.circular==1,'circular must be boolean or 0-1')


    % MBB (fix) block length
    assert(numel(opt.lblock)==1 && opt.lblock>1,'lblock must be a number >1 [N of hours]')

    dt = nanmax(data.t) - nanmin(data.t);  % need the statistics toolbox
    % dt = max(data.t(~isnan(data.t))) - min(data.t(~isnan(data.t))); % [days]
    if opt.lblock > 0.1*(dt*24) % the block length must be at least 1/10 of series length in hours
        warning('check_boot_opt: block length greater than 1/10 of the observed period')
    end

    % MBB (variable) block length distribution name [block length sampler]:
    legal_distrib = {'unif' 'geom' 'fix'}; % 'pois'
    opt.lblock_dist = lower(opt.lblock_dist);
    assert(any(strcmp(opt.lblock_dist,legal_distrib)), ['Not a valid lblock_dist for MBB, possible distributions are ' legal_distrib{:}])


    % MBB (variable) block length distribution parameter(s)
    if isempty(opt.lblock_par) 

        % set default block length parameters for the chosen block dist
        switch opt.lblock_dist
            case 'geom'
                opt.lblock_par = 1/(30*24 +1); % probability parameter p
            case 'unif'
                opt.lblock_par = 60*24;   % max block length in hours 
                % opt.lblock_par = [0, 60*24];   % min-max block length in hours
            case 'fix'
                opt.lblock_par = opt.lblock;
            % OLD: removed since --> Gaussian with this parameter magnitude    
            % case 'pois' 
            %     opt.lblock_par = 30*24 + 1;
        end
    else

        % check the block length parameters for the chosen block dist
        switch opt.lblock_dist
            case 'geom'
                assert(numel(opt.lblock_par )==1,'MBB block length parameter (geometric distribution) must be a scalar (1/mean block length)')
                assert(opt.lblock_par>0 && opt.lblock_par<1,'MBB block length parameter must be a probability in (0,1)')
            case 'unif' 
                assert(numel(opt.lblock_par)==2,'MBB block length vector must be a 2 element vector (min-max block lengths for an uniform distribution)')
            % case 'pois' 
            %     assert(numel(opt.lblock_par)==1,'MBB block length parameter (poisson distribution) must be a scalar (mean block length)') 
            case 'fix'
                assert(numel(opt.lblock_par)==1,'MBB block length must be a scalar')
                assert(~isnan(opt.lblock_par), 'block length cannot be NaN for fix length MBB')
                assert(opt.lblock_par>1 && opt.lblock_par<opt.lblock_par,'MBB Block length must be in (1,T]')
        end

    end
    
   
    assert(isempty(opt.iboot),'cannot provide a vector of time indices to be resampled for method=mbb')
    
    % remove the SPB options and other unusefull opt
    opt = rmfield(opt,'noise');
    if ~strcmpi(opt.lblock_dist,'fix')
        opt = rmfield(opt,'lblock','iboot'); 
    end

end


function opt = check_spb_opt(opt)

    legal_psd   = {'fft','fftnoise','skfft', 'sknoise', 'skewed','cpfft', 'cpnoise','dllftt', 'dllnoise', 'daniell'};
    legal_names = {'fft','fftnoise','skfft', 'sknoise', 'cpfft', 'cpnoise','dllftt', 'dllnoise'};

    % check the name of the methods to simulate the residuals
    opt.noise = lower(opt.noise);
    assert(any(strcmp(opt.noise,legal_psd)), ['Not a valid noise value for SPD, possible methods are ' legal_names{:}])
    if ~any(strcmp(opt.noise,{'fft','fftnoise'}))
        warning([opt.noise ' SPB presently deprecated > use the ftt option instead'])
    end

    % % remove MMB options (the block info)
    assert(isempty(opt.iboot),'cannot provide a vector of time indices to be resampled for method=spb')
    opt = rmfield(opt,{'lblock','lblock_par','lblock_dist','circular','block_output','iboot','iblocks'});

end


function opt = check_knownit_opt(opt)

        % check that a vecrot of time indices is provided
        assert(~isempty(opt.iboot),'a vector of time indices must be provided for empty method')
        
        % check the provided blocks and time indices
        lt = length(data.t);
        nb = opt.nboot;
        [r,c] = size(opt.iboot);
        assert(r==lt, ['iboot must be a ' num2str(lt) ' x ' num2str(nb) ' matrix. Provided is ' num2str(r) ' x ' num2str(c) '.'])
        if c~=nb
            warning(['Decrasing nboot to ' num2str(c) ' to match the provided ibott size'])
        end
        
        assert(all(isinteger(opt.iboot(:))), 'iboot must be a matrix of integers (time vector indices)')
        assert( all(opt.iboot(:)>1 & opt.iboot(:)<lt), [ 'iboot must contain indices in [1,' num2str(lt) ']' ])
         
        
        % % remove MMB and SPB options 
        opt = rmfield(opt,{'lblock','lblock_par','lblock_dist','block_output','iblocks'});
        
end