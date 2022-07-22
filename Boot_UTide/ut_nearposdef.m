
%%--------------------------------------------------------- 
function covposdef = ut_nearposdef(cov,maxit)
    % UT_NEARPOSDEF()
    % find nearest positive definite matrix by Higham 2002 algorithm
    %   with identity weight matrix
    % inputs
    %   cov - matrix for which to find the nearest positive definite matrix
    %   maxit - max iterations before quitting (optional; default 1000)
    % output
    %   covposdef - the nearest positive definite covariance matrix
    %             - diag(diag(cov)) if max iterations limit is met 
    % UTide v1p0 9/2011 d.codiga@gso.uri.edu
    % (builds on validcorr.m by Erland Ringstad, by adding:
    %   conversion from correlation matrix to covariance matrix;
    %   end of search once all eigenvalues are positive;
    %   optional maxit input. )
    
    if size(cov,1)~=size(cov,2)
        fprintf('ut_nearposdef: cov must be square and real.\n');
        covposdef = nan;
        return;
    end
    
    if all(eig(cov)>=0)
        covposdef = cov;
        return;
    end
    
    % set the max num of iterations
    if ~exist('maxit','var')
        maxit = 1000;    
    elseif maxit<=0 || rem(maxit,1)
           fprintf('ut_nearposdef: maxit must be positive integer.\n');
           covposdef = nan;
        return;
    end
    
    S = zeros(size(cov));
    Y = cov;
    k = 0;
    while  ~all(eig(Y)>=0) && k<maxit
            R = Y - S;
        [Q,D] = eig(R);
            X = Q*max(D,0)*Q';
            S = X - R;
            Y = X - diag(diag(X)-1);
            k = k + 1;
    end
    
    if ~all(eig(Y)>=0)
        %warning('ut_nearposdef:diagSig','ut_nearposdef: reached max iterations (%d); zeroing off-diagonals.')
        % fprintf(sprintf(['ut_nearposdef: reached max iterations (%d); zeroing off-diagonals.\n'],maxit));
        covposdef = diag(diag(cov));
    else
              std = sqrt(diag(cov));
        covposdef = Y.*(std*std');
    end
    end
    