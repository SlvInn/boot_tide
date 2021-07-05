%%--------------------------------------------------------- 
function bopt = cut_bootinit(varargin)
% S. Innocenti. 2019/09 - modified 2021/06

% set defaults for the bootstrap parameters
bopt.seed      = 2*floor(rand(2,1)*(10^8)) +1; % define an odd number as pseudo-numebr generator seed 
bopt.kw        = coef.W; % use HA wls/irls weights as knwon weights

bopt.mtd       = 'MBB';  % bootstrap type
bopt.circular  = false;  % use circular MBB or not
bopt.lBlk_dist = 'geom'; % default distribution (name)
% % % OLD:
% % % bopt.lBlk      = 31;     % mean block length [days] = 56 days for daily data
% % % bopt.lBlk_dist = 'pois'; % default distribution (name)
% % % bopt.lBlk      = 1/31;   % mean block length [days] = 31 days for daily data
% % % NEW: bopt.lBlk  defined according to the selected dist (see below)
bopt.lBlk      = [] 

% read user-defined optional arguments
optargin = size(varargin,2);

i = 1;
while i <= optargin
    switch lower(varargin{i})
        case {'method' 'mtd'}
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
           assert(numel(bopt.seed)==2,'bootstrap random seed must be a 2 element vector')
           
        case {'knownweights' 'kweights' 'kw'}
           bopt.kw = varargin{i+1};
           assert(numel(bopt.kw )==ltin || isempty(bopt.kw ),'known weight vector must be empty or have the same length as uin,vin')
           
        otherwise 
            error(['cut_boot1: unrecognized input: ' varargin{i}]); 
           
    end
    i = i + 2;
end

% set the default block length or parameter, if not specified by the user
if isempty(bopt.lBlk) 
    switch bopt.lBlk_dist
        case 'geom'
            bopt.lBlk= 1/31
        case 'unif' 
            bopt.lBlk= [0,60]
        case 'pois' 
            bopt.lBlk= 30
        case 'fix'
            bopt.lBlk= 30
    end
end    


% controls on user-defined inputs
% assert(~(isempty(bopt.lBlk) &&  strcmpi(bopt.mtd,'mbb')), 'you must specify min-max block length (''blocklength'') for the MBB method with non-hourly data ') 

