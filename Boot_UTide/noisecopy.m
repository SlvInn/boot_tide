function noise = noisecopy(x,Nseries)
% Generate noise with a given power spectrum while also controlling 
% skewnees and the kurtosis of the simulated signal
%
% INPUTS:
% x: the observed time series (a column vector)
% Nseries: number of noise series to generate. (default=1)
% 
% OUTPUT:
% noise: surrogate series with same power spectrum as f. (each column is a surrogate). 
% --->> Silvia Innocenti, adapted from Aslak Grinsted (2009)
if nargin<2
    Nseries=1;
end


% FFT
x  = x(:);
f  = fft(x);
N  = length(x);


% number of the unique fft points
Np = floor((N-1)/2);
% disp('Np')
% disp(Np)

% generate a random phase
% phases = rand(Np,Nseries)*2*pi ;
% 
% 
% % phases in the freq domain signal
% phik = complex(cos(phases),sin(phases)); % this was the fastest alternative in my tests. 
% % phik2 = exp(j*phases ); 

% test: gaussian complex phase
% s2x = var(x(:));
% phik = sqrt(s2x/2)*(randn(Np,Nseries)+j*randn(Np,Nseries));

% generate a random phase
phases = randn(Np,Nseries);
% phases = sqrt(1/2)*(randn(Np,Nseries)+j*randn(Np,Nseries));
m    = 0;
s    = 1; %sqrt(nanvar(x));
sk   = skewness(x(~isnan(x)));
ku   = kurtosis(x(~isnan(x)));
phases = m + s.*phases+ sk.*(phases.^2)+ku.*(phases.^3);

phik =  -pi + (2*pi)*((phases-nanmin(phases))/(nanmax(phases)-nanmin(phases)));
phik = exp(j*phik);

% use the same fft for all simulations
f  = repmat(f,1,Nseries);

% frequency domain signal(s)             
% f(2:Np+1,:)          = f(2:Np+1,:).*phik  ;
% f(end:-1:end-Np+1,:) = conj(f(2:Np+1,:));
% 1, 2:Np+1, 51, N:-1:N-Np+1
if rem(N, 2) %if even
    zk = f(2:Np,:).*phik(1:end-1,:);
    zk = [f(1,:); zk; f(Np+1,:);  f(Np+2,:); flipud(conj(zk));];
else  %if odd
    zk = f(2:Np+1,:).*phik;
    zk = [f(1,:); zk; f(Np+2,:); flipud(conj(zk));];
end

% disp('rFFTnoise')
% disp(size(zk))
% disp([1,Np+1;  N-Np+1,N])

% produce the noise
% noise =real(ifft(f,[],1)); 
noise = real(ifft(zk)); 

end
