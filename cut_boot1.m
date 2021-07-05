%%--------------------------------------------------------- 
function [coef_boot,B] = cut_boot1(tin,uin,vin,coef,n_boot,varargin)
% S. Innocenti 
% 2020/04/25
%
% Bootstrap for a single record. 
% See comments and help for cut_boot(), cut_renconstr(), and cut_solv().
%
% TO DO (warning messages have been set up for most of these issues):
% - finish controls on user-defined inputs
% - test the fucntions for 2D analyses; 
% - 

% heritate estimation options
opt = coef.aux.opt;

% reconstrauct series 
[uhat,vhat] = cut_reconstr1(tin,coef,'cnstit',coef.name);

% check series dimension
if opt.twodim   
    nDim = 2;
    Yhat = [uhat,vhat];
    Yres = [uhat-uin,vhat-vin];
    warning('the code must be checked and tested for the two dimentional case')   
else 
    nDim = 1;
    Yhat = uhat;
    Yres = uhat-uin;   
end


% get the input length and remore NaNs
ltin = numel(tin);
if opt.twodim
    uvgd = ~isnan(uin) & ~isnan(vin);
else
    uvgd = ~isnan(uin);
end
t = tin(uvgd);
nt = numel(t);
     
% define the LOR for equispaced or irregularly sampled series        
if opt.equi == 1 % based on times; u/v can still have nans ("gappy")
    lor = (max(tin)-min(tin)); 
else
    lor = (max(t) - min(t));   
end
    

% Bootstrap options >> Set defaults
bopt.mtd       = 'MBB';  % bootstrap type
bopt.circular  = false;  % use circular MBB or not
bopt.lBlk      = 31;     % mean block length [days] = 56 days for daily data
bopt.lBlk_dist = 'pois'; % default distribution (name)
bopt.seed      = 2*floor(rand(2,1)*(10^8)) +1; % define an odd number as pseudo-numebr generator seed 
bopt.kw        = coef.W; % use HA wls/irls weights as knwon weights


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

% controls on user-defined inputs
% assert(~(isempty(bopt.lBlk) &&  strcmpi(bopt.mtd,'mbb')), 'you must specify min-max block length (''blocklength'') for the MBB method with non-hourly data ') 

%% Prepare matrix for regression, based on hat-resonstructions
         
% define empty inference
% % assert(isempty(opt.infer),'Inference is not yer considered in the bootstrap function')

% Get the basis function
% options to compute the basis function
ngflgs = [coef.aux.opt.nodsatlint coef.aux.opt.nodsatnone coef.aux.opt.gwchlint coef.aux.opt.gwchnone];

% compute basis function 
[E,F,U,V] = cut_E(t,coef.aux.reftime,coef.aux.frq,coef.aux.lind,coef.aux.lat,ngflgs,coef.aux.opt);
        
if ~opt.complexf % if sin/cos formulation 
    B = [F.*cos(2*pi*(U+V)), F.*sin(2*pi*(U+V))];
    
    if ~isempty(opt.infer)
        error 'not done yet: inference must be adapted to the sin/cos definition of the model >> try with complexf option = true'
        % see the other part of the IF statement (when opt.infer is ~empty)
    else
          Q = [];
       beta = [];
         nR = 0;
        nNR = numel(coef.name);
    end
else % if complex formulation 
    B = [E conj(E)];
    
    if ~isempty(opt.infer)
        error 'not done yet: inference must be adapted to the sin/cos definition of the model >> try with complexf option = true'
    else
          Q = [];
       beta = [];
         nR = 0;
        nNR = numel(coef.name);
    end
    
end

% include (a) column(s) for temporal trends
if opt.notrend
    B = [B ones(nt,1)];
else
    B = [B ones(nt,1) (t-coef.aux.reftime)/lor];   
end
    
% define a weighted matrix, if needed
if ~isempty(bopt.kw)
    w = sqrt(bopt.kw(uvgd));
    W = diag(w);
    wB = (W*B); % OK
else
    w = [];
    wB = B;
end         
          
           
 % empties 
 coef_boot.name = coef.name;    
          nM    = numel(coef.M);
 coef_boot.M    = nan(nM,n_boot);
          nC    = numel(coef.name);
 coef_boot.g    = nan(nC,n_boot);
if ~opt.complexf                     
    coef_boot.A    = nan(nC,n_boot);
    coef_boot.mean = nan(1,n_boot);
    if ~opt.notrend
        coef_boot.slope = nan(1,n_boot);
    end        
else
    coef_boot.Lsmaj = nan(nC,n_boot);
    coef_boot.Lsmin = nan(nC,n_boot);
    coef_boot.theta = nan(nC,n_boot);
    coef_boot.umean = nan(1,n_boot);
    coef_boot.vmean = nan(1,n_boot);
    if ~opt.notrend
        coef_boot.uslope = nan(1,n_boot);
        coef_boot.vslope = nan(1,n_boot);
    end
end
 

% internal variables
   dt0e = 24*(tin(end) - tin(1)); % number of hourly time steps
 expNbl = ceil(4*(dt0e/(24*bopt.lBlk))); 

 
 tic 
if strcmpi(bopt.mtd,'mbb')
    
    % % construct bootstrap resamples 
    
    % empties
   lBlocks = cell(n_boot,1);
   nBlocks =  nan(n_boot,1);
    I0boot = cell(n_boot,1);
    ITboot = cell(n_boot,1);
    
    switch bopt.lBlk_dist
        
        case 'unif'
    
        assert(numel(bopt.lBlk)==2,'MBB Block length vector must be a 2 element vector (min-max block lengths for an Uniform distribution)')
        
        % min-max block length in hours
        minlBlk = bopt.lBlk(1)*24;
        maxlBlk = bopt.lBlk(2)*24;
        
        rng(bopt.seed(2))
        lBlkrand = randi([minlBlk maxlBlk],expNbl,n_boot);
    

        case 'pois' 
        
        assert(numel(bopt.lBlk)==1,'MBB Block length parameter (poisson distribution) must be a scalar (mean block lengths)')
        
        meanlBlk = bopt.lBlk;
            rng(bopt.seed(2))
        lBlkrand = 24.*poissrnd(meanlBlk-1, expNbl,n_boot)+1;

        case 'geom' 
        
            assert(numel(bopt.lBlk)==1,'MBB Block length parameter (geometric distribution) must be a scalar (1/mean block lengths)')
            assert(bopt.lBlk>0 && bopt.lBlk<1,'MBB Block length parameter must be a prob in (0,1)')

            meanlBlk = bopt.lBlk;
            rng(bopt.seed(2))
        lBlkrand = 24.*geornd(meanlBlk, expNbl,n_boot);

        case  'fix'

        assert(numel(bopt.lBlk)==1,'MBB Block length (fixed length) must be a scalar')
        assert(bopt.lBlk>=1 && bopt.lBlk<floor(ltin/24),'MBB Block length must be in [1,T]')
         
        meanlBlk = bopt.lBlk;
        lBlkrand = 24.*repmat(meanlBlk, expNbl,n_boot);
   end
   
   % extract random starts of the blocks
    rng(bopt.seed(1))
    I0rand = randi([1 ltin], expNbl,n_boot);
    
   
   
%    if strcmpi(bopt.lBlk_dist,'fix') && bopt.lBlk>=floor(ltin/24)-1
%         I0rand = 1;
%    else
%     rng(bopt.seed(1))
%     I0rand = randi([1 ltin], expNbl,n_boot);
%    end
  
   if bopt.circular
       tcirc = [tin(:); tin(:)];
       Ycirc = [Yres;Yres];
   else
       tcirc = tin(:);
       Ycirc = Yres;
   end
    
   % resamples
     for bs = 1 : n_boot
         
        Tbs = [];
       iTbs = [];
     lBlkbs = [];
         bl = 0;
         
        % print the iteration steps 
        %if mod(bs,10^3)==0
        %   disp([' ' num2str(bs) 'th boot iteration - ' datestr(now)])
        %end


        while numel(Tbs)< ltin
            bl = bl+1;
                   
            i0 = I0rand(bl,bs);  % random start of block
        lBlk = lBlkrand(bl,bs);

             tE = (tin(i0).*24+lBlk)./24;
             iE = find(tcirc<=tE,1,'last');
             iE = max([iE,i0]); 

            Tbs = [Tbs; tcirc(i0:iE)];
         lBlkbs = [lBlkbs; numel(i0:iE)];  
           iTbs = [iTbs; [i0:iE]'];

        end

               
               iTbs = iTbs(1:ltin);
        lBlocks{bs} = lBlkbs;
        nBlocks(bs) = bl;
         I0boot{bs} = I0rand(1:bl);
         ITboot{bs} = iTbs;

        % resamples
        Yboot = nan(ltin,nDim);
        for d = 1 : nDim
            Yboot(:,d) = Yhat(:,d) + Ycirc(iTbs,d);%Yres(iTbs,d);
        end
        
        
            yboot = Yboot(:,1);
        if  opt.twodim
            yboot = complex(Yboot(:,1),Yboot(:,2));
        end

            
        % estimate model 
        m_boot = boot_reg(wB, yboot(uvgd), w, opt);
        coef_boot.M(:,bs) = m_boot;
            
     end
    

      % update the bootstrap structure 
     coef_boot.boot.iT0boot  = I0boot;
     coef_boot.boot.ITboot   = ITboot;
     coef_boot.boot.lBlocks  = lBlocks;
     coef_boot.boot.nBlocks  = nBlocks;
     
else
    
    % remove block infos from the bootstrap structure
    bopt.lBlk = [];        
    bopt.lBlk_dist = [];
    tmmpseed = bopt.seed;
    
    % empty matrix of seeds
    pbootseed = nan(n_boot,nDim);

            % interpolate missing values
            FTTres = nan(nt,nDim);
            
            if opt.equi == 1  
                
                % eliminate Nan
                % nonna  = find(~any(isnan(Yres),2)); 
                % lfft   = length(nonna);
                % FTTres = nan(nt,nDim);
                nonna = find(uvgd);
                
                for d = 1 : nDim
                    
                    %   % interpolate missing values - OLD METHOD
                    %            na = find( isnan(Yres(:,d)));   
                    %           nna = find(~isnan(Yres(:,d))); 
                    %    Yres(na,d) = interp1(nna,Yres(nna,d),na,'pchip',0);
                    %    % print(mean(Yres))


                    % estimate residual PSD 
                    FTTres(:,d) = fft(Yres(nonna,d)); 
                end
            else
                
                warning('For irregularly spaced observations, residuals are interpolated over a time vector with fixed time step')
                teq = linspace(min(t), max(t), nt);
                
                for d = 1 : nDim
                    % interpolate non-missing values over a regular
                    % time vector [to be changed]
                    YresEq = interp1(t,Yres(uvgd,d),teq,'pchip',0);
                
                    % estimate residual PSD
                    FTTres(:,d) = fft(YresEq);  
                end

            end

            % geenerate random seeds for each boot    
            rng(tmmpseed(d))
            pbootseed(:,d) = 2*floor(rand(n_boot,1)*(10^5)) +1;
            
            % save residua PSD
            coef_boot.boot.ResPsd = nan(nt,1);
            coef_boot.boot.ResPsd(uvgd,:) = FTTres;
            
            
            % % construct bootstrap resamples 
             for bs = 1 : n_boot
                 
                % print the iteration steps
%                  if mod(bs,10^3)==0
%                     disp([' ' num2str(bs) 'th boot iteration - ' datestr(now)])
%                  end

                   % construct the resamples 
                    Yboot = nan(ltin,nDim);
                    for d = 1 : nDim
                        
                        rng(pbootseed(bs,d)) % set the seed
                        BootErr  = fftnoise(FTTres(:,d),1); % simulate a noise
                        Yboot(uvgd,d) = Yhat(uvgd,d) + BootErr ; % sum the simulated noise to the original reconstruction
                    end
         
                    % define the yboot sample at this iteration
                    yboot = Yboot(:,1);
                    if  opt.twodim
                       yboot = complex(Yboot(:,1),Yboot(:,2));
                    end
                       

                     % estimate model 
                      m_boot = boot_reg(wB, yboot(uvgd), w, opt);
           coef_boot.M(:,bs) = m_boot;                           
           
                    
                    
             end 
            
end
% from M to A,g, mean, and trend:
coef_boot = cut_from_mboot_to_param(coef_boot,n_boot,opt,nNR,nR);
coef_boot.boot.comp_time = toc;


% output
coef_boot.boot_opt = bopt;
coef_boot.boot_opt.aux = coef.aux;
end


            


function m_boot = boot_reg(wB, yboot, w, opt) 
% Apply the least square estimation on the model defined by the B basis
% function and with the given opt 

    switch opt.method

        case 'ols'
             iva = ~isnan(yboot);
             m_boot = wB(iva,:)\yboot(iva);

        case {'cauchy','andrews','bisquare','fair','huber', 'logistic','talwar','welsch','ridge', 'lasso', 'enet'}

            if any(strcmpi(opt.method,{'ridge', 'lasso', 'enet'}))
                warning([opt.method ' bootstrap not implemented yet >> using known weights'])
                assert (~isempty(w), [opt.method ' bootstrap need known weights'] )
            end

            if ~isempty(w)  
                % use the knwon weights in the regression  
                   iva = ~isnan(yboot) & ~isnan(w);
                m_boot = wB(iva,:)\(w(iva).*yboot(iva));
            else
                % re-estimate weights
                   iva = ~isnan(yboot);
                m_boot = robustfit(wB(iva,:),ctranspose(yboot(iva)),opt.method,opt.tunconst,'off');
            end

     end   
end

