function coef_boot = cut_bootcore(coef,tin,uvgd,B,Yres,Yhat,bopt,n_boot)

opt = coef.aux.opt;
if opt.twodim   
    nDim = 2;
else 
    nDim = 1;
end

% define a weighted matrix, if needed
if ~isempty(bopt.kw)
    w  = sqrt(bopt.kw(uvgd));
    wB = (diag(w)*B); % OK
else
    w  = [];
    wB = B;
end  

%% --  define empties -- %%
coef_boot.name = coef.name;    
nM   = numel(coef.M);
coef_boot.M    = nan(nM,n_boot);
nC   = numel(coef.name);
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
ltin = numel(tin); % time length
dt0e = 24*(tin(end) - tin(1)); % number of hourly time steps
expNbl = ceil(4*(dt0e/(24*bopt.lBlk))); 

% time vector without nan:
t  = tin(uvgd);
nt = length(t)

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

coef_boot.boot.comp_time = toc;

end



