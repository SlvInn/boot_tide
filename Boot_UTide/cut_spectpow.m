function [coef,Puu,Puv,Pvv] = cut_spectpow(coef,opt,e,t,tin,tgd,uvgd,elor)
% spectral power (Puu, Pvv, Puv) of residual
% fprintf('conf. int''vls ... ');

   % band-averaged (ba) spectral densities
    if opt.equi % if the time vector ....
        if sum(tgd)>sum(uvgd) % if the num of ~nan(time) > the num of ~nan uv   
              efill = interp1(t,e,tin(tgd));
			  
            if any(isnan(efill)) % fill start&/end nans w/ nearest good
                        ind = find(isnan(efill));
                       ind2 = ind(ind<find(~isnan(efill),1,'first')); 
                efill(ind2) = efill(max(ind2)+1);
                       ind2 = ind(ind>find(~isnan(efill),1,'last'));
                efill(ind2) = efill(min(ind2)-1);
            end
			
            ba = ut_pdgm(tin(tgd),efill,coef.aux.frq,1,0); % band-averages of line-decimated periodograms         
        else
            ba = ut_pdgm(tin(tgd),e,coef.aux.frq,1,0); % band-averages of line-decimated periodograms
        end
    else
            ba = ut_pdgm(t,e,coef.aux.frq,0,opt.lsfrqosmp); % band-averages of line-decimated periodograms
    end
   
	
    % power [ (e units)^2 ] from spectral density [ (e units)^2 / cph ]
            df = 1/(elor*24); % [cph/elor]
        ba.Puu = ba.Puu*df;   % [units^2/cph * df] = [units^2/elor] ?!
    if opt.twodim
        ba.Pvv = ba.Pvv*df;
        ba.Puv = ba.Puv*df;
    end
    
    % assign band-avg power values to NR & R freqs
        Puu = zeros(size(coef.aux.frq));
    if opt.twodim
        Pvv = Puu;
        Puv = Pvv;
    else
        Pvv = [];
        Puv = [];
    end
    
    for i = 1 : length(ba.Puu)
                 ind = find(coef.aux.frq>=ba.fbnd(i,1) & coef.aux.frq<=ba.fbnd(i,2));
            Puu(ind) = ba.Puu(i);
			
        if opt.twodim
            Pvv(ind) = ba.Pvv(i);
            Puv(ind) = ba.Puv(i);
        end
    end
	if ~opt.twodim
		coef.aux.rpds = Puu; % output Residual PSD. Added (SI) 
	else
	    coef.aux.rpds{1} = Pvv; % output Residual PSD. Added (SI)
		coef.aux.rpds{2} = Puv; % output Residual PSD. Added (SI)
    end
end



%%--------------------------------------------------------- 
function P = ut_pdgm(t,e,cfrq,equi,frqosmp)
% UT_PDGM()
% spectral densities (band-averages of line-decimated periodograms)
% inputs
%          t = times [datenum UTC] (nt x 1)
%          e = error (residual; model-raw misfit), complex/real 2D/1D (nt x 1)
%       cfrq = frequencies of NR & R constituents [cycles per hour] (nc x 1)
%       equi = 0/1 if equispaced times or not (1 x 1)
%   frqosamp = lomb-scargle freq oversampling factor (ignored if equi=1)
% outputs (all 9 x 1)
%                   P.Puu = 1-sided auto-sp dens of u err (=real(e)) [e units^2/cph] (9 x 1)
%                  P.fbnd = edges of freq bands [cph] (9 x 2) 
%   P.Pvv (2dim case only)= as P.Puu but for v err (=imag(e)) (9 x 1) 
%   P.Puv (2dim case only)= cross-sp dens between u and v (9 x 1)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (band-averaging and line decimating from t_tide, Pawlowicz et al 2002)

    nt = length(e);
    hn = hanning(nt);
    if equi
        [Puu1s,allfrq] = pwelch(real(e),hn,0,nt);
    else
        [Puu1s,allfrq] = ut_lmbscga(real(e),t,hn,frqosmp);
    end
               fac = (nt-1)/(2*pi*(t(end)-t(1))*24); % conv fac: rad/sample to cph
            allfrq = allfrq*fac; % to [cycle/hour] from [rad/samp]
             Puu1s = Puu1s/fac; % to [e units^2/cph] from [e units^2/(rad/samp)]
    [P.Puu,P.fbnd] = ut_fbndavg(Puu1s,allfrq,cfrq);  
     % disp(P)
     
    if ~isreal(e)
        if equi
            [Pvv1s,~] = pwelch(imag(e),hn,0,nt);
            [Puv1s,~] = cpsd(real(e),imag(e),hn,0,nt);
        else
            [Pvv1s,~] = ut_lmbscga(imag(e),t,hn,frqosmp);
            [Puv1s,~] = ut_lmbscgc(real(e),imag(e),t,hn,frqosmp);
        end

            Pvv1s = Pvv1s/fac;
        [P.Pvv,~] = ut_fbndavg(Pvv1s,allfrq,cfrq);
            Puv1s = real(Puv1s)/fac;
        [P.Puv,~] = ut_fbndavg(Puv1s,allfrq,cfrq);
            P.Puv = abs(P.Puv);

    end
end

%%--------------------------------------------------------- 
function [Pxx,F] = ut_lmbscga(x,t,w,ofac)
% UT_LMBSCA()
% Auto-spectral density estimate based on the mean-removed Lomb-Scargle 
%   periodogram; designed such that, in the special case of equispaced
%   input times and no frequency oversampling, the output units,
%   normalization, and frequency values replicate those of pwelch(). 
%
% Inputs: 
%   * Sequence of real values x from n arbitrarily distributed times t
%   * Taper/window shape w at n equispaced times (e.g. w=hanning(n) )
%   * Frequency oversampling factor ofac (integer >=1)
%
% Outputs:
%   * One-sided auto-spectral density estimate Pxx [units: (x-units^2) /
%      (radians per sample)] as column vector, based on mean-removed and
%      unnormalized Lomb-Scargle periodogram, after applying window/taper w
%      to entire record (segmenting is not implemented); at zero and
%      Nyquist frequencies Pxx is computed by the FFT-based auto-spectral
%      density formula.
%   * Frequencies F (units: radians per sample) of the Pxx estimates, as 
%      column vector: for ofac=1, these are the same equispaced values as
%      for an FFT-based spectral calculation using a time series of n
%      values at uniformly distributed times from min(t) to max(t), thus F
%      has resolution 2*pi/n (rad/sample) and for even n it spans the
%      interval [0,Pi] with (n+2)/2 values, while for odd n it spans the
%      interval [0,Pi) with (n+1)/2 values; for ofac>1, the same ranges are
%      spanned but with resolution increased by ofac times.
%
% Notes:
%   Designed such that, for the special case of input values x from 
%      uniformly distributed times t, [Pxx,F]=ut_lmbscga(x,t,w,1) 
%      gives results (both Pxx and F) identical to [Pxx,F]=pwelch(x,w,0,n).
%   The cross-periodogram companion function is ut_lmbscgc().
%   Finally, for context with respect to the now-deprecated psd() 
%      function, the functions pwelch() and psd() are related as follows.
%      Consider [Pxxw,Fw] = pwelch(x,w,0,n) and [Pxxp,Fp] = psd(x,n,1/dt),
%      for real x, w = hanning(n), and sampling interval dt [t units]. 
%      A. Pxxw is the one-sided power spectral density with units 
%      [x-units^2/(rad/sample)] and Pxxp is the two-sided power spectrum
%      (not spectral density) with units [x-units^2]. Therefore, to obtain
%      Pxxw from Pxxp requires (i) division by the sampling frequency
%      (1/dt) because Pxxp is power spectrum not power spectral density,
%      yielding dt*Pxxp [units: x-units^2/(cycles per t-unit)]; (ii)
%      division by 2*pi and by dt, to convert from cycles per t-unit to
%      radians per sample, and (iii) multiplication by 2 to convert from
%      two-sided to one-sided, for frequencies other than 0 and Nyquist. 
%       Even n: Pxxw is [Pxxp(1)/(2*pi) Pxxp(2:end-1)/pi Pxxp(end)/(2*pi)];
%       Odd n:  Pxxw is [Pxxp(1)/(2*pi) Pxxp(2:end)/pi].
%      B. Units of Fw are radians per sample and units of Fp are cycles 
%      per t-unit, so to obtain Fw from Fp requires multiplying by 2*pi 
%      and by dt. Even or odd n: Fw is 2*pi*dt*Fp.
% UTide v1p0 9/2011 d.codiga@gso.uri.edu  

    % reshape and check inputs
    x = x(:)'; 
    t = t(:)';
    w = w(:)';
    if ~isequal(size(x),size(t),size(w)) || ~(isreal(x)&&isreal(t)&&isreal(w)) 
        error('ut_lmbscga: x,t, and w must be same size and real.');
    end
    ofac = round(ofac);
    if ofac<1
        error('ut_lmbscga: ofac must be >= 1.');
    end
    n = length(x);
    % taper
    dt = (max(t) - min(t))/(n-1);
     w = interp1(min(t):dt:max(t),w,t);
    xw = x.*w;
    % frequencies (excluding 0 and nyquist)
    F = ((1/ofac):(1/ofac):n/2)'*2*pi/n; % radian/sample
    if ~mod(n,2) % even, exclude nyqyuist
        F(end) = [];
    end
    % unnormalized, 2-sided, mean-removed lomb-scargle periodogram [x-units^2] 
    arg = (2/dt)*(F*t)';
    tau = 0.5*imag(log(sum(cos(arg))+1i*sum(sin(arg)))); % see angle() hdr
    arg = 0.5*arg - tau(ones(n,1),:);
     xb = mean(xw);
    Pxx = sum(sparse(1:n,1:n,xw-xb)*cos(arg)).^2./sum(cos(arg).^2);
    Pxx = 0.5*(Pxx+sum(sparse(1:n,1:n,xw-xb)*sin(arg)).^2./sum(sin(arg).^2));
    % include zero-frequency
    Pxx = [n*xb^2; Pxx']; 
      F = [0; F];
    if ~mod(n,2)  % even, include nyquist
          c = [ones(1,n/2); -1*ones(1,n/2)];
        Pxx = [Pxx; n*mean(xw.*c(:)')^2];
          F = [F; pi];
    end
    % recover taper suppression
    Wss = n*sum(w.*w);
    Pxx = n^2*Pxx/Wss;

    % divide by 2*pi to get spectral density rel radians/sample
    Pxx = Pxx/(2*pi); % two-sided, [x-units^2/(rad/sample)]
    % multiply by 2 to get one-sided (except zero and nyquist)

    Pxx(2:end-1) = 2*Pxx(2:end-1); % one-sided, [x-units^2/(rad/sample)]
    if mod(n,2) % odd n, highest frequency not nyquist
        Pxx(end) = 2*Pxx(end); % one-sided, [x-units^2/(rad/sample)]
    end
end

%%--------------------------------------------------------- 
function [Pxy,F] = ut_lmbscgc(x,y,t,w,ofac)
% UT_LMBSCGC()
% Cross-spectral density estimate based on the mean-removed Lomb-Scargle 
%   cross-periodogram; designed such that, in special case of equispaced
%   input times and no frequency oversampling, the output units,
%   normalization, and frequency values replicate those of CPSD().  
%
% Inputs: 
%   * Real-valued sequences x and y from n arbitrarily distributed times t
%   * Taper/window shape w at n equispaced times (e.g. w=hanning(n) )
%   * Frequency oversampling factor ofac (integer >=1)
%
% Outputs:
%   * One-sided cross-spectral density estimate Pxy [units: (x-units)
%      (y-units)/(radians per sample)], based on the mean-removed and
%      unnormalized Lomb-Scargle cross-periodogram, after applying
%      window/taper w to the entire records (segmenting is not
%      implemented); at zero and Nyquist frequencies Pxy is computed by the
%      FFT-based cross-spectral density formula.
%   * Frequencies F (units: radians per sample) of the Pxy estimates: for 
%      ofac=1, these are the same equispaced values as for an FFT-based 
%      spectral calculation using a time series of n values at uniformly
%      distributed times from min(t) to max(t), thus F has resolution
%      2*pi/n (rad/sample) and for even n it spans the interval [0,Pi] with
%      (n+2)/2 values, while for odd n it spans the interval [0,Pi) with
%      (n+1)/2 values; for ofac>1, the same ranges are spanned but with
%      resolution increased by ofac times.
%
% Notes:
%   Designed such that, for the special case of input values x and y from 
%      uniformly distributed times t, [Pxy,F]=ut_lmbscgc(x,y,t,w,1) 
%      gives results (both Pxy and F) identical to [Pxx,F]=cpsd(x,y,w,0,n).
%   The auto-periodogram companion function is ut_lmbscga(). 
%   The relationship between cpsd() and the deprecated csd() function is
%      fully analogous to the relationship between pwelch() and the
%      deprecated psd() function; an explanation of the latter is given at
%      the end of the function description text in ut_lmbscga.m.
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

    % reshape and check inputs
    x = x(:)';
    y = y(:)';
    t = t(:)';
    w = w(:)';
    if ~isequal(size(x),size(y),size(t),size(w)) || ...
            ~(isreal(x)&&isreal(y)&&isreal(t)&&isreal(w)) 
        error('ut_lmbscgc: x,y,t, and w must be same size and real.');
    end
    ofac = round(ofac);
    if ofac<1
        error('ut_lmbscgc: ofac must be >= 1.');
    end
    n = length(x);
    % taper
    dt = (max(t) - min(t))/(n-1);
     w = interp1(min(t):dt:max(t),w,t);
    xw = x.*w;
    yw = y.*w;
    % frequencies (excluding 0 and nyquist)
    F = ((1/ofac):(1/ofac):n/2)'*2*pi/n; % radian/sample
    if ~mod(n,2) % even, exclude nyqyuist
        F(end) = [];
    end
    % unnormalized, 2-sided, mean-removed lomb-scargle cross-periodogram 
    % [(x-units)(y-units)] 
    arg = (1/dt)*(F*t)';
    tau = 0.5*imag(log( sum(cos(2*arg)) + 1i*sum(sin(2*arg)) ));
     xb = mean(xw);
     yb = mean(yw);
     X1 = sum(sparse(1:n,1:n,xw-xb)*cos(arg));
     X2 = sum(sparse(1:n,1:n,xw-xb)*sin(arg));
     Y1 = sum(sparse(1:n,1:n,yw-yb)*cos(arg));
     Y2 = sum(sparse(1:n,1:n,yw-yb)*sin(arg));
    arg = arg - tau(ones(n,1),:);
      A = 1./sqrt(sum(cos(arg).^2));
      B = 1./sqrt(sum(sin(arg).^2));
     F0 = exp(1i*tau);
      X = F0.*(A.*X1 + 1i*B.*X2);
      Y = conj(F0).*(A.*Y1 - 1i*B.*Y2);
    Pxy = 0.5*X.*Y;
    % include zero-frequency
    Pxy = [n*xb*yb; Pxy'];
    F = [0; F];
    if ~mod(n,2)  % even, include nyquist freq
        c = [ones(1,n/2); -1*ones(1,n/2)];
        Pxy = [Pxy; n*mean(xw.*c(:)')*mean(yw.*c(:)')];
        F = [F; pi];
    end
    % recover taper suppression
    Wss = n*sum(w.*w);
    Pxy = n^2*Pxy/Wss;

    % divide by 2*pi to get spectral density rel radians/sample
    Pxy = Pxy/(2*pi); % two-sided, [x-units^2/(rad/sample)]
    % multiply by 2 to get one-sided (except zero and nyquist)
    Pxy(2:end-1) = 2*Pxy(2:end-1); % one-sided, [x-units^2/(rad/sample)]
    if mod(n,2) % odd n, highest frequency not nyquist
        Pxy(end) = 2*Pxy(end); % one-sided, [x-units^2/(rad/sample)]
    end
end

%%--------------------------------------------------------- 
function [avP,fbnd] = ut_fbndavg(P,allfrq,cfrq)
% UT_FBNDAVG()
% line-decimate and band-average spectra
% inputs
%        P = periodogram to treat [e units^2/cph]
%   allfrq = frequency values of (equispaced) P estimates [cph]
%     cfrq = frequencies of constituents [cph] (nc x 1)
% outputs
%    avP = line-decimated and band-averaged spectrum [e units^2/cph] (9 x 1)
%   fbnd = frequencies [cph] at edges of averaged bands (9 x 2)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu
% (based on residual_spectrum.m of t_tide, Pawlowicz et al 2002)

  df = allfrq(3)-allfrq(2);
P(round(cfrq./df)+1) = NaN; 
fbnd = [.00010 .00417;
        .03192 .04859;
        .07218 .08884;
        .11243 .12910;
        .15269 .16936;
        .19295 .20961;
        .23320 .25100;
        .26000 .29000;
        .30000 .50000];
nfbnd = size(fbnd,1);
  avP = zeros(nfbnd,1);
for k = nfbnd:-1:1
    b1 = find(allfrq>=fbnd(k,1));
    b2 = find(allfrq<=fbnd(k,2));
    b3 = find(isfinite(P));
    jbnd = intersect(intersect(b1,b2),b3); 
    if any(jbnd)
        avP(k) = mean(P(jbnd));
    elseif k<nfbnd
        avP(k) = P(k+1);   
    end
end
end
