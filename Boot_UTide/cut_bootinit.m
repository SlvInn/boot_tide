%%--------------------------------------------------------- 
function bopt = cut_bootinit(coef,tin,varargin)
% S. Innocenti. 2019/09 - modified 2021/06
% Set the default values of the options and retrieve user-defined parameters.
%
% {OPTIONS/VARARGIN}:
% 'method'/'mtd', mtdname - string indicating which bootstrap algorithm must be used 
%               to construct the residual resamples. Accepted strings are:
%          'mbb' for implementing the Moving Block Boostrap (either with fixed or random block length')
%          'parboot', 'pboot', or 'psdboot' for the Semi-Parametric Bootstrap based on the residual power spectrum
%           Default: 'mbb'
%
% 'circular', true/false  - boolean flag determining if a circular strategy 
%           [Politis and Romano (1994)] must be used in constructing the resamples.
%           Default: false
%
% 'lblk_dist'/'blockdist', distname - string defining the probability distribution 
%            to extract/simulate the block lengths for each MBB resample. 
%            Allowed distname:
%            * 'fix'  - blocks with fixed length = ndays (see the 'blocklength' input below), 
%                       no distribution involved. 
%            * 'geom' - blocks with random lengths esimulated from a geometric distribution with 
%                       parameter 1/p = ndays (see the 'blocklength' input below)
%            * 'pois' - blocks with random lengths esimulated from a poisson distribution with 
%                       parameter lambda = ndays (see the 'blocklength' input below) 
%            * 'unif' - blocks with random lengths esimulated from an uniform distribution with 
%                       parameters [a,b] = ndays (see the 'blocklength' input below)
%             Default: blockdist = 'geom'
%          
% 'blocklength'/'lblk', ndays - block length [days] for fixed length block in MBB 
%                               Default: ndays = 31
%
% 'lblkpar'/'blocklengthparam', theta - block length distribution parameter(s) when using 
%           random blocks MBB. Specifically, theta has the following interpretations:
%            * inverse of the average block length for blockdist ='geom'. Default: 1/31 
%            * average block length for blockdist = 'pois'. Default: 30 
%            * min and max block length for blockdist = 'unif' (i.e. ndays is a 2-element vector). Default: [0,60] 
%
%
% NOTE: The 'circular', 'blockdist', 'blocklength', and 'blocklengthparam' options
%       are ignored when using SPB. 
%
%
% 'seed', seeds - 2-element vector containing the random seeds to be used in the 
%           resampling. For MBB, the first seed is used in the simulation of the block 
%           random starts, while the second is used in the simulation of the block lengths.
%           For the SPB, only the first seed is used to initialize the fttnoise function.  
%
% 'knownweights'/'kweights'/'kw', kwvect - Tx1 or empty vector of known WLS weights that will be used  
%       in the HA regression of each resample. Default: kwvect = coef.W 
%       (i.e. use the original HA wls/irls weights, if present in coef). 
%       IF kwvect = [], the IRLS method is applied for each resample (robustfit.m), 
%       otherwise, the known weights are used as in a GLS. 
%       These weights represent estimates of the reciprocal of the observation std 
%       NOTE: weights are assumed w = 1 / sigma_ii, while some IRLS models use
%       bisquared weights [i.e., w = 1 / (sigma_ii^2)]; consider
%       using the sqrt of the weights through the option "knownweights", if needed.       
%       
% 'knownresiduals'/'yres'/'residuals', kyres - Tx1 vector or Tx2 matrix  of known regression residuals 
%       to be used to construct the bootstrap resamples. Default: [] - Estimated by cut_reconstr1.m
%
% 'knownreconstr'/'yhat'/'reconstr' , kyhat - Tx1 vector or Tx2 matrix  of known regression predictions  
%       to be used to construct the bootstrap resamples. Default: [] - Estimated by cut_reconstr1.m
%
%

% set default variables for resampling the regression model
bopt.kw        = coef.W;  % use the HA wls/irls weights as knwon weights
bopt.kyhat      = [];     % tidal reconstruction  based on the estimated HA model
bopt.kyres      = [];     % HA model residuals 
% bopt.kB        = [];    % known basis function for the model

% set defaults for the bootstrap parameters
bopt.seed      = 2*floor(rand(2,1)*(10^8)) +1; % define an odd number as pseudo-numebr generator seed 
bopt.mtd       = 'MBB';  % bootstrap type
bopt.circular  = false;  % use circular MBB or not
bopt.lBlk_dist = 'geom'; % default distribution (name)
bopt.lBlk      = [] ; % Block length or Block length distribution parameter >> bopt.lBlk is defined according to the selected dist (see below)

% read user-defined optional arguments
optargin = size(varargin,2);

i = 1;
while i <= optargin
   
    switch lower(varargin{i})
        case {'method','mtd'}
            bopt.mtd = varargin{i+1};
            assert(any(strcmpi(bopt.mtd,{'mbb' 'pboot' 'parboot' 'psdboot'})),'not valid bootstrap method')
          
        case {'circular'}
            bopt.circular = varargin{i+1};
            assert(numel(bopt.circular)==1,'circular must be scalar (boolean or 0-1)')
            assert(islogical(bopt.circular) || bopt.circular==0 || bopt.circular==1,'circular must be boolean or 0-1')
                                
        case {'lblk' 'blocklength' 'lblkpar' 'blocklengthparam'}
            bopt.lBlk = varargin{i+1};
               
        case {'lblk_dist' 'blockdist'}
            bopt.lBlk_dist = varargin{i+1};
            legal_distrib = {'unif' 'pois' 'geom' 'fix'};
            assert(any(strcmpi(bopt.lBlk_dist,legal_distrib)),['Not a valid sitribution for the simulation of MBB Block lengths,possible distributions are ' legal_distrib{:}])
       
        case 'seed' 
           bopt.seed = varargin{i+1};
           assert(numel(bopt.seed)==2,'the bootstrap random seed must be a 2 element vector')
           
        case {'knownweights' 'kweights' 'kw'}
           bopt.kw = varargin{i+1};
           % OLD: 
           % assert(numel(bopt.kw )==ltin || isempty(bopt.kw ),'known weight vector must be empty or have the same length as uin,vin')
           % NEW:
           if ~isempty(bopt.kw)
                assert(numel(bopt.kw)==numel(coef.W),'the input weight vector must have the same length as the estimated vector')
           end

        % case {'knownbasis' 'basis' 'b'}   
        %     bopt.kB = varargin{i+1};
        %     if ~isempty(bopt.kB)
        %         [rB,cB] = size(B)
        %         assert(rB==numel(coef.W),'basis matrix size is not consistent with observation and time vector length')
        %         assert(cB==numel(coef.M),'basis matrix size is not consistent with the number of estimated coefficients in coef.M')
        %     end

        case {'knownresiduals' 'yres' 'residuals'}
            bopt.kyres  = varargin{i+1}; 
            if coef.aux.opt.twodim 
                assert(size(bopt.kyres,2)==2,'residuals must be a 2-column matrix')
                assert(size(bopt.kyres,1)==numel(tin),'residuals and time vector sizes are not consistence')    
            else    
                assert(length(bopt.kyres)==numel(tin), 'residuals must have the same length as observation and time vectors')
            end

        case {'knownreconstr' 'yhat' 'reconstr'}    
            bopt.kyhat  = varargin{i+1};  
            
            if coef.aux.opt.twodim 
                assert(size(bopt.kyhat,2)==2,'HA reconstruction must be a 2-column matrix')
                assert(size(bopt.kyhat,1)==numel(tin),'HA reconstruction and time vector sizes are not consistence')    
            else    
                assert(length(bopt.kyhat)==numel(tin), 'reconstructed series must have the same length as observation and time vectors')
            end


        otherwise 
            error(['cut_boot1: unrecognized input: ' varargin{i}]); 
           
    end
    i = i + 2;
end

% set the default block length or parameter, if not specified by the user
if isempty(bopt.lBlk) 
    switch bopt.lBlk_dist
        case 'geom'
            bopt.lBlk= 1/31;
        case 'unif' 
            bopt.lBlk= [0,60];
        case 'pois' 
            bopt.lBlk= 30;
        case 'fix'
            bopt.lBlk= 30;
    end
end    


% controls on user-defined inputs
% assert(~(isempty(bopt.lBlk) &&  strcmpi(bopt.mtd,'mbb')), 'you must specify min-max block length (''blocklength'') for the MBB method with non-hourly data ') 

if ~isempty(bopt.kyres)
    assert( ~isempty(bopt.kyhat),'HA reconstruction are needed in input when using user-defined residuals')
end
if ~isempty(bopt.kyhat)
    assert( ~isempty(bopt.kyres),'known residuals are needed when providing user-defined HA reconstructions')
end

end