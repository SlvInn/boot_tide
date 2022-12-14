function out = get_spb(data,opt)
% Get the simulated y resamples with the Semi-Parametric Bootstrap (SPB) method.
%
% INPUT:
% data  - structure of processed data series   >> See help check_data
% opt   - structure with the bootstrap options >> See help boot_tide.m for a description of the opt structure fields 
% 
% OUTPUT -  structure with (some of) the following fields:
% res_ftt    - (T x S) matrix of the estimated residual spectra
% y_boot     - (T x nboot x S) matrix of the nboot simulated resamples for each spatial location S
% theta_boot - (P x nboot x S)  matrix of the nboot replicates of tide model parameters estimated on the resamples
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2022/08

    % time and residual vector without nan
    t    = data.t(data.valid);
    yres = data.yres(data.valid,:);


    % number of spatial locations
    ns   = size(yres,2); 

    % n of resamples
    nb = opt.nboot;

    % geenerate and store random seeds for each boot 
    rng(opt.seed(1))
    seeds = 2*floor(rand(nb,ns)*(10^5)) +1;

    % To handle irregularly spaced observations, residuals are interpolated over a time vector with fixed time step
    lt  = numel(data.t);
    teq = roundn(linspace(min(t), max(t), lt),-4); % rounded at the hour resolution

    % define the size of the result structure 
    if isfield(opt,'tide_model')
        nout = numel(opt.theta); 
    else
        nout = lt;      
    end

    % empties
    out.res_ftt = nan(lt,ns); % empty matrix for the residual FTT 
    res = nan(nout,nb,ns); % empty matrix of y boot resamples or boot param
    for s = 1 : ns % for each location 

        % interpolate non-missing values over a regular time vector [to be changed]
        yreseq = interp1(t,yres(:,s),teq,'pchip',0); %'pchip'

        % generate nb random residuals based on the fft and store residual FFT
        [e_boot_s, out.res_ftt(:,s)] = get_spb_noise(yreseq , opt.noise, nb);
        % TODO: verify the possibility of simulating all the series for nb >= 10^4 
        % >> if not possible, break the get_spb_noise functions (using the fft in input)

        % % for each resample, estimate the tidal parameters using opt.tide_model or get the yboot reconstr
        res(:,:,s) = from_spb_res_to_out(data.yhat(:,s),e_boot_s, opt);
        
    end


    % put the results in the relevant field of the out structure
    if isfield(opt,'tide_model')
        out.theta_boot = res;
    else
        out.y_boot = res;
    end    

end

function out = from_spb_res_to_out(y,es,opt) 

    nb = size(es,2);
    if isfield(opt,'tide_model') % compute the parameters 
       
        % % for each resample, estimate the tidal parameters using opt.tide_model
        for b = 1 : nb
            yb = y + es(:,b);
            out(:,b) = opt.tide_model(yb);
        end
    else

        % % get all the resamples of the original series
        for b = 1 : nb
            out(:,b) = y + es(:,b);
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               NOISE GENERATORS                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions that generate noise series with same spectrum 
% as an observed series x

function [noise,f] = get_spb_noise(y_residual,method,n_boot)
% generate the spb residual resamples using one noise generator
% INPUTS:
% y_residual: time series (column vector).
%     method: string indicating the noise generator (funtion) name,
%             either fft, fftnoise, skfft, sknoise, cpfft, cpnoise, dllftt, or dllnoise. 
% n_boot: number of noise series to generate. 
% 
% OUTPUT:
%  noise: simulated noise series with same power spectrum as f (each column is a resample).
%      f: spectrum estimate used by the noise generator
    
    switch method
        case {'fft','fftnoise'}
            [noise,f] = fftnoise(y_residual,n_boot);
        case {'skfft', 'sknoise', 'skewed'}    
            [noise,f] = sknoise(y_residual,n_boot);
        case {'cpfft', 'cpnoise'}
            [noise,f] = cpnoise(y_residual,n_boot);
        case {'dllftt', 'dllnoise', 'daniell'}   
            [noise,f] = dllnoise(y_residual, n_boot);
    end

end




function [noise,f] = fftnoise(x,n_series)
% Generate noise with the same power spectrum as the observed x series.
%  
% noise=fftnoise(f[,n_series])
%
% INPUTS:
% x: time series (column vector)
% n_series: number of noise series to generate. (default=1)
% 
% OUTPUT:
%  noise: surrogate series with same power spectrum as f (each column is a simulation).
%      f: fast fourier transform of the observed series.
%
% --- Adapted from Aslak Grinsted (2009)
% https://www.mathworks.com/matlabcentral/fileexchange/32111-fftnoise-generate-noise-with-a-specified-power-spectrum?s_tid=srchtitle
%

    if nargin<2
        n_series=1;
    end

    % f  = f(:); % AG - Original
    x  = x(:);   % SI
    f  = fft(x); % SI
    N  = length(f);
    Np = floor((N-1)/2);
    phases = rand(Np,n_series)*2*pi;

    phases = complex(cos(phases),sin(phases)); % AG: this was the fastest alternative in my tests. 



              ff = repmat(f,1,n_series);
    ff(2:Np+1,:) = ff(2:Np+1,:).*phases;
    ff(end:-1:end-Np+1,:) = conj(ff(2:Np+1,:));

    noise = real(ifft(ff,[],1)); 

end


function  [noise,ff] = dllnoise(x,n_series)
% Generate noise with a given power spectrum.
% Useful helper function for Monte Carlo null-hypothesis tests and confidence interval estimation.
%  
% noise=fftnoise(x,n_series])
%
% INPUTS:
%       x: time series (column vector)
% n_series: number of noise series to generate. (default=1)
% 
% OUTPUT:
%  noise: surrogate series with same power spectrum as f (each column is a surrogate).
%     ff: Daniell avg the spectrum estimated via the fast fourier transform of the observed series.
%
% --->> Silvia Innocenti, adapted from Aslak Grinsted (2009)

    if nargin<2
        n_series=1;
    end

    x  = x(:);  
    f  = fft(x); 
    N  = length(f);
    Np = floor((N-1)/2);


    % Daniell avg the spectrum (boxcar avg)
    B   = ones(16,1) ;
    B   = B./sum(B) ;
    msk = ones(N,1) ;

    
    % apply a running mean of 16 values
    ff = conv2(f,B,'same')./conv2(msk,B,'same'); %for older Matlab versions
    

    % simulate the random phases
    phases = rand(Np,n_series)*2*pi;
    phases = complex(cos(phases),sin(phases)); % this was the fastest alternative in my tests. 

    % multiply the phases and the ftt for each resample
                    f = repmat(ff,1,n_series);
            f(2:Np+1,:) = f(2:Np+1,:).*phases;
            
    % apply the IFFT         
    f(end:-1:end-Np+1,:) = conj(f(2:Np+1,:));
    noise = real(ifft(f,[],1)); 

end

function [noise,ff] = sknoise(x,n_series)
% Generate noise with a given power spectrum while also controlling 
% skewnees and the kurtosis of the simulated signal.
%
% INPUTS:
% x: the observed time series (a column vector)
% n_series: number of noise series to generate. (default=1)
% 
% OUTPUT:
%  noise: surrogate series with same power spectrum as f (each column is a simulation).
%     ff: fast fourier transform of the observed series.
%
% --->> Silvia Innocenti, adapted from Aslak Grinsted (2009)

    if nargin<2
        n_series=1;
    end


    % FFT
    x  = x(:);
    ff = fft(x);
    N  = length(x);


    % number of the unique fft points
    Np = floor((N-1)/2);



    % phases in the freq domain signal
    % test: gaussian complex phase
    % s2x = var(x(:));
    % phik = sqrt(s2x/2)*(randn(Np,n_series)+j*randn(Np,n_series));
    % phik = complex(cos(phases),sin(phases)); %


    % generate a random phase
    phases = randn(Np,n_series);
    % phases = sqrt(1/2)*(randn(Np,n_series)+j*randn(Np,n_series));

    m    = 0;
    s    = 1; 
    % s = sqrt(nanvar(x)); % test
    % s = nanvar(x); % test
    sk   = skewness(x(~isnan(x)));
    ku   = kurtosis(x(~isnan(x)));
    phases = m + s.*phases+ sk.*(phases.^2)+ku.*(phases.^3);

    phik =  -pi + (2*pi)*((phases-nanmin(phases))/(nanmax(phases)-nanmin(phases)));
    phik = exp(j*phik);

    % use the same fft for all simulations
    f  = repmat(ff,1,n_series);


    % frequency domain signal(s)             
    if rem(N, 2) %if even
        zk = f(2:Np,:).*phik(1:end-1,:);
        zk = [f(1,:); zk; f(Np+1,:);  f(Np+2,:); flipud(conj(zk));];
    else  %if odd
        zk = f(2:Np+1,:).*phik;
        zk = [f(1,:); zk; f(Np+2,:); flipud(conj(zk));];
    end


    % produce the noise
    noise = real(ifft(zk)); 

end




function [noise,f] = cpnoise(series,n_series)
% Copy noise: generate noise with the same (red) spectrum as the observed x series.
%
% INPUTS:
% x: the observed time series (a column vector)
% n_series: number of noise series to generate. (default=1)
% 
% OUTPUT:
%  noise: surrogate series with same power spectrum as f (each column is a simulation).
%      f: fast fourier transform of the observed series.
% 
% --->> Silvia Innocenti, adapted from rednoise.m - H. Zhivomirov  (2013)

    series = series(:);
    m = length(series);

    % define the length of the noise vector and ensure  
    % that M is even, this will simplify the processing
    N = round(m); 
    if rem(N, 2)
        M = N+1;
    else
        M = N;
    end

    % prepare a vector with frequency indexes 
    NumUniquePts = M/2 + 1;  % number of the unique fft points
    k = 1:NumUniquePts;     % vector with frequency indexes 

    % interpolate missing valies
            na = find(isnan(series));   
           nna = find(~isnan(series));   
    series(na) = interp1(nna,series(nna),na,'pchip');
          f = fft(series);
                    

    noise = nan(M,n_series);
    for s = 1 : n_series

        % generate a white noise sequence
        x = randn(1, M);

        % FFT
        X = fft(x);


        % manipulate the left half of the spectrum so the PSD
        % is proportional to the PSD of the observed series
        X = X(1:NumUniquePts);  
        Xobs = f(1:NumUniquePts)';  
        X = X.*Xobs;

        % prepare the right half of the spectrum - a conjugate copy of the left
        % one except the DC component and the Nyquist component - they are unique,
        % and reconstruct the whole spectrum
        X = [X conj(X(end-1:-1:2))];

        % IFFT
        y = real(ifft(X));

        % ensure that the length of y is N
        y = y(1, 1:N);

        % form the noise matrix and ensure unity standard 
        % deviation and zero mean value (columnwise)
        y = reshape(y, [m, 1]);
        y = bsxfun(@minus, y, mean(y));
        noise(:,s) = bsxfun(@rdivide, y, std(y));

    end    

end