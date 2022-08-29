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
    assert(numel(opt.lblock)==1 && opt.lblock>1,'lblock must be a number >1 [N of days]')

    dt = nanmax(data.t) - nanmin(data.t);
    if opt.lblock > 0.1*dt
        warning('check_boot_opt: block length greater than 1/10 of the observed time period')
    end
else
    opt.lblock = NaN;
end    

# MBB (variable) block length distribution name [block length sampler]:
legal_distrib = {'unif' 'pois' 'geom' 'fix'};
assert(any(strcmpi(opt.lblock_dist,legal_distrib)), ['Not a valid lblock_dist for MBB,possible distributions are ' legal_distrib{:}])

# MBB (variable) block length distribution parameter(s)
if strcmpi(opt.method,'mbb')
   if isempty(opt.lblock_par) 
        disp(['set default block length parameter for ' opt.lblock_dist 'MBB'])
        switch opt.lBlk_dist
            case 'geom'
                opt.lblock_par = 1/31;
            case 'unif' 
                opt.lblock_par = [0,60];
            case 'pois' 
                opt.lblock_par = 30;
            case 'fix'
            assert(~isnan(opt.lblock),"block length cannot be NaN for fix length MBB")
            opt.lblock_par = opt.lblock
        end
        
   end
else
    opt.lblock_par = NaN;
end

% SILVIA: SONO QUI >> finisci di fare cut_bootinit.m from line 166
% aggiungi dei commenti e dei warnings per i parametri messi a Nan 
end