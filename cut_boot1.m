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

% reconstruct the original series based on the estiamted parameters (i.e. regression fit)
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
    

% Bootstrap options >> Set defaults and read the ser-defined inputs
bopt = cut_bootinit(varargin)

%% Prepare matrix for regression, based on hat-resonstructions
         
% define (empty) inference (TO DO)
% % assert(isempty(opt.infer),'Inference is not yer considered in the bootstrap function')

% Get the basis function
% options to compute the basis function
ngflgs = [coef.aux.opt.nodsatlint coef.aux.opt.nodsatnone coef.aux.opt.gwchlint coef.aux.opt.gwchnone];

% compute the model basis function 
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

 
tic % start the time counteer
if strcmpi(bopt.mtd,'mbb')
    
    % % construct bootstrap resamples 
    
    % empties
   lBlocks = cell(n_boot,1);
   nBlocks =  nan(n_boot,1);
    I0boot = cell(n_boot,1);
    ITboot = cell(n_boot,1);

    % simulate the block lengths (if random bloks):
    lBlkrand = cut_boot_blk_lengt(bopt)
    
   
   % simulate/draw random starts of the blocks
    rng(bopt.seed(1))
    I0rand = randi([1 ltin], expNbl,n_boot);
    
   
% % (SI): OLD >>   
% %    if strcmpi(bopt.lBlk_dist,'fix') && bopt.lBlk>=floor(ltin/24)-1
% %         I0rand = 1;
% %    else
% %     rng(bopt.seed(1))
% %     I0rand = randi([1 ltin], expNbl,n_boot);
% %    end
  
   if bopt.circular
       tcirc = [tin(:); tin(:)];
       Ycirc = [Yres;Yres];
   else
       tcirc = tin(:);
       Ycirc = Yres;
   end
    
   % construct resamples
     for bs = 1 : n_boot

        % construct a series of obs indices to select data in the resamples
         
        % empties for the bs-th resample
        Tbs = []; % series  of the resampled dates  
       iTbs = []; % indices of the resampled dates  
     lBlkbs = []; % lengths of each resampled block
         bl = 0;  % block counter

        while numel(Tbs)< ltin 
            bl = bl+1; % increase the block counter
                   
            i0 = I0rand(bl,bs);  % random start of the bl-th block
        lBlk = lBlkrand(bl,bs);  % length of the bl-th block

             tE = (tin(i0).*24+lBlk)./24;   % end time of the bl-th block
             iE = find(tcirc<=tE,1,'last'); % index of the end time in the series
             iE = max([iE,i0]);  

            Tbs = [Tbs; tcirc(i0:iE)];   % attach the block of dates to Tbs
         lBlkbs = [lBlkbs; numel(i0:iE)];% record the length of the present block  
           iTbs = [iTbs; [i0:iE]'];      % attach the block of indices to iTbs

        end

        % construct the resamples by summing the resamples residual 
        % and the original reconstruction
        iTbs = iTbs(1:ltin);
        Yboot = nan(ltin,nDim);
        for d = 1 : nDim
            Yboot(:,d) = Yhat(:,d) + Ycirc(iTbs,d);
        end
        
            yboot = Yboot(:,1);
        if  opt.twodim
            yboot = complex(Yboot(:,1),Yboot(:,2));
        end

        % store the simulated quantities
        lBlocks{bs} = lBlkbs;
        nBlocks(bs) = bl;
         I0boot{bs} = I0rand(1:bl);
         ITboot{bs} = iTbs;
            
        % estimate the regression parameters for the resample 
        m_boot = cut_boot_reg(wB, yboot(uvgd), w, opt);
        coef_boot.M(:,bs) = m_boot; % store the paraemeters in coef_boot
            
     end
    

      % update the bootstrap structure 
     coef_boot.boot.iT0boot  = I0boot;
     coef_boot.boot.ITboot   = ITboot;
     coef_boot.boot.lBlocks  = lBlocks;
     coef_boot.boot.nBlocks  = nBlocks;
     
else
    
    % remove the block info from the bootstrap structure
    bopt.lBlk = [];        
    bopt.lBlk_dist = [];
    tmmpseed = bopt.seed;
    
    % empty matrix of seeds
    pbootseed = nan(n_boot,nDim);

            % empty matrix for the residual FTT 
            FTTres = nan(nt,nDim);
            
            if opt.equi == 1  
                nonna = find(uvgd); %locate NaNs
                for d = 1 : nDim 
                    % estimate the residual PSD at non-NaN time steps
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
            
            
            % % construct the parametric bootstrap resamples 
             for bs = 1 : n_boot

                % construct the resamples 
                Yboot = nan(ltin,nDim);
                for d = 1 : nDim
                    
                    rng(pbootseed(bs,d)) % set the seed
                    BootErr  = fftnoise(FTTres(:,d),1); % simulate a noise with same spectrum as the residuals
                    Yboot(uvgd,d) = Yhat(uvgd,d) + BootErr ; % sum the simulated noise to the original reconstruction
                end
        
                % define the yboot sample at this iteration
                yboot = Yboot(:,1);
                if  opt.twodim
                    yboot = complex(Yboot(:,1),Yboot(:,2));
                end
                    

                % estimate the regression parameters for the resample 
                m_boot = cut_boot_reg(wB, yboot(uvgd), w, opt);
     coef_boot.M(:,bs) = m_boot;                           
             
             end 
            
end

% from M to A,g, mean, and trend:
coef_boot = cut_from_mboot_to_param(coef_boot,n_boot,opt,nNR,nR);
coef_boot.boot.comp_time = toc;

% output
coef_boot.boot_opt = bopt; % bootstrap options
coef_boot.boot_opt.aux = coef.aux; % HA model setup
end


            



