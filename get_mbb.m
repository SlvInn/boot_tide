function out = get_mbb(data,opt)
% Get the sequential indices to be used for constructing the y resamples
% with the Moving Block Bootstrap (MBB) method.
%
% INPUT:
% data  - structure of processed data series   >> See help check_data
% opt   - structure with the bootstrap options >> See help boot_tide.m for a description of the opt structure fields 
% 
% OUTPUT -  structure with (some of) the following fields:
% i_boot     - (T x nboot) matrix of time indices to use for constructing the nboot resamples of y [i.e., s.t. y_boot = yhay + yres(i_boot,:,:)]
% blocks     - (nboot x 2) cell with the series of block lengths and block start indices for each bootstrap resample
% y_boot     - (T x nboot x S) matrix of the nboot simulated resamples for each spatial location S
% theta_boot - (P x nboot x S)  matrix of the nboot replicates of tide model parameters estimated on the resamples
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2022/08

    % time and residual vector without nan
    t    = data.t(data.valid);
    lt   = numel(data.t);
    yres = data.yres(data.valid,:);
    ns = size(yres,2);

    % n of resamples
    nb = opt.nboot; 


    % define the size of the result structure 
    if isfield(opt,'tide_model')
        nout = numel(opt.theta); 
    else
        nout = lt;  
    end  

    % empties
    i_boot  = nan(lt,nb);  % index 
    blocks = cell(nb,2);  % index of the start time and length of each block of each bootstrap resample
    res = nan(nout,nb,ns); % empty matrix of y boot resamples or boot param

    % simulate the block lengths (if random blocks, otherwise repeat the fix length) and the random starts of the blocks:
    [rand_l_blk, rand_i0, time_h] = mbb_block_lengths(t,data.tresol,opt);

    % double the time vector to allow for circular bootstrap
    if opt.circular
        tc = [t(:); t(:)];
    else
        tc = t(:);
    end 
    
    % construct the series of obs indices to select data in the resamples
    for b = 1 : nb

        % empties for the b-th resample
        % tb  = []; % series  of resampled times/dates
        itb   = []; % indices of the resampled times/dates
        lblk  = []; % lengths of each block
        blk   = 0 ; % block counter

        while numel(itb) < lt
           blk = blk+1; % increase the block counter
            i0 = rand_i0(blk,b);     % random start of the blk-th block
            l  = rand_l_blk(blk,b);  % length of the blk-th block

            tE = (time_h(i0)+l)/24;       % end time of the blk-th block [in days]
            iE = find(tc<=tE,1,'last'); % index of the end time in the original series
            iE = max([iE,i0]);  


          % tb  = [tb; tc(i0:iE)];    % attach the block of dates to tb
          lblk  = [lblk; numel(i0:iE)];  % record the length of the block  
          itb   = [itb; [i0:iE]'];       % attach the block of indices to itb

        end

        % eliminate exceeding indices and store the results in the output variables
        i_boot(:,b) = itb(1:lt);
        blocks{b,1} = lblk(:);
        blocks{b,2} = rand_i0(1:blk,b);
        res(:,b,:)  = from_mbb_res_to_out(data,i_boot(:,b),opt);
    end


    % add the block construction infor to the output structure, if requested
    if opt.block_output
        out.i_boot = i_boot;
        out.blocks = blocks;
    end    

    % put the results in the relevant field of the out structure
    if isfield(opt,'tide_model')
        out.theta_boot = res;
    else
        out.y_boot = res;
    end

end 


function [rl_blk, rand_i0, th] = mbb_block_lengths(t,tresol,opt)
% For each MBB block, simulate the block length or get a vector of (equal) fixed lengths

    % length of the vector of valid time steps
    nt    = length(t);

    nh   = 24./tresol;        % n of hours 
    mt_h = round(min(t)*nh); % min time in hours: [days]*nh
    Mt_h = round(max(t)*nh); % max time in hours: [days]*nh
    th = (mt_h:1:Mt_h)';     % hourly time vector [h]
    lth = numel(th);         % N of hourly time steps
    
    % number of possible block length values (more than needed)
    if strcmp(opt.lblock_dist,'geom')
     n_blocks = ceil(10*lth*opt.lblock_par);
    else
        n_blocks = ceil(10*lth./opt.lblock_par); % 10 x expected number of blocks 
    end

    
    rng(opt.seed(2)) % use the second random seed
    switch opt.lblock_dist
        case 'geom'
            rl_blk = geornd(opt.lblock_par, n_blocks, opt.nboot);

        case 'unif' 
            rl_blk = randi([1 opt.lblock_par], n_blocks, opt.nboot);
          
        case 'pois' 
            rl_blk = poissrnd(opt.lblock_par-1, n_blocks, opt.nboot)+1;
          
        case 'fix'
            rl_blk = repmat(opt.lblock_par, n_blocks, opt.nboot);
    end
 
    % simulate/draw random starts of the blocks
    rng(opt.seed(1))
    rand_i0 = randi(nt, n_blocks, opt.nboot);
end




function out = from_mbb_res_to_out(data,ib,opt) 
% TODO: check if this funtion works with varous locations because the tidal model shoudl be the same at all locations

    % reconstruction and resampled residuals
    y  = data.yhat;
    ns = size(y,2);
    eb = data.yres(ib,:);


    if isfield(opt,'tide_model') % compute the parameters 
        
        
        % % for each resample, estimate the tidal parameters using opt.tide_model
        for s = 1 : ns
            yb = y(:,s) + eb(:,s);
            out(:,s) = tide_model(yb)
        end
    else

        % % get all the resamples of the original series
        out = y + eb;
        
    end
end
