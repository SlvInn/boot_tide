function [a,pha] = from_bsincos_to_amppha(m)
    % convert sin-cos amplitudes to tidal amplitudes and phases
    
    % n of tidal constituents
     m = m(:);
    nc = (length(m)-1)/2; % assume to have the intercept in m
    
    ap = m(1:nc);
    am = m(nc+(1:nc));
    
    a   = sqrt(ap.^2 + am.^2);
    pha = atan2d(am,ap);
    
    
    % ut_cluster
          I = ones(1,size(pha,2));
        ii  = (pha-pha(:,I ))> 360/2;
    pha(ii) = pha(ii)-360;
         ii = (pha-pha(:,I))< - 360/2;
    pha(ii) = pha(ii)+360;
    
    pha(pha<0) = 360+pha(pha<0);
    
    end