function opt = check_boot_opt(data,opt)
% Apply basic checks on the entered bootstrap options
%
% INPUT:
# data - structure of processed data series   >> See help check_data
% opt  - structure with the bootstrap options >> See help boot_tide.m for a description of the opt structure fields 
% 
% OUTPUT:
% opt  - structure with the bootstrap options >> See help boot_tide.m
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2020/04 - 2022/08

# number of bootstrap replicates/resamples
assert(isinteger(opt.nboot),'nboot must be an integer')

# random seed for resamplers
assert(numel(bopt.seed)==2,'the bootstrap random seed must be a 2 element vector')

# bootstrap method: Moving Block Bootstrap (MBB) or Semi-Parametric Bootstrap (SPB)
assert(any(strcmpi(opt.method,{'mbb' 'spb'})),'not valid bootstrap method')

# MBB circular option (avoid truncation at series end)
assert(numel(opt.circular)==1,'circular must be scalar (boolean or 0-1)')
assert(islogical(opt.circular) || opt.circular==0 || opt.circular==1,'circular must be boolean or 0-1')
   
# MBB (fix) block length
if strcmpi(opt.method,'mbb') 
    assert(numel(opt.lblock)==1 && opt.lblock>1,'lblock must be a number >1 [N of hours]')

    dt = nanmax(data.t) - nanmin(data.t); % [days]
    if opt.lblock > 0.1*(dt*24) % the block length must be at least 1/10 of series length in hours
        warning('check_boot_opt: block length greater than 1/10 of the observed time period')
    end
else

    % remove the block info from the bootstrap structure
    opt.lblock      =  NaN;
    opt.lblock_par  =  NaN;
    opt.lblock_dist = 'none';
    
end    

# MBB (variable) block length distribution name [block length sampler]:
legal_distrib = {'unif' 'pois' 'geom' 'fix'};
assert(any(strcmpi(opt.lblock_dist,legal_distrib)), ['Not a valid lblock_dist for MBB,possible distributions are ' legal_distrib{:}])

# MBB (variable) block length distribution parameter(s)
if strcmpi(opt.method,'mbb')
    if isempty(opt.lblock_par) 
        ## disp(['set default block length parameter for ' opt.lblock_dist 'MBB'])

        switch opt.lBlk_dist
            case 'geom'
                opt.lblock_par = 1/(30*24 +1); % probability parameter p
            case 'unif' 
                opt.lblock_par = [0, 60*24];   % min-max block length in hours
            case 'pois' 
                opt.lblock_par = 30*24 + 1;
            case 'fix'
                opt.lblock_par = opt.lblock
        end
    else

        switch opt.lBlk_dist
            case 'geom'
                assert(numel(opt.lblock_par )==1,'MBB block length parameter (geometric distribution) must be a scalar (1/mean block length)')
                assert(opt.lblock_par>0 && opt.lblock_par<1,'MBB block length parameter must be a probability in (0,1)')
            case 'unif' 
                assert(numel(opt.lblock_par)==2,'MBB block length vector must be a 2 element vector (min-max block lengths for an uniform distribution)')
            case 'pois' 
                assert(numel(opt.lblock_par)==1,'MBB block length parameter (poisson distribution) must be a scalar (mean block length)') 
            case 'fix'
                assert(numel(opt.lblock_par)==1,'MBB block length must be a scalar')
                assert(~isnan(opt.lblock_par),"block length cannot be NaN for fix length MBB")
                assert(opt.lblock_par>1 && opt.lblock_par< opt.lblock_par,'MBB Block length must be in (1,T]')
        end




    end
else
    opt.lblock_par = NaN;
end

end