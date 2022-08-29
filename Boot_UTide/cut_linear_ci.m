function [coef,varcov_mCw, varcov_mCc] = cut_linear_ci(opt,varcov_mCw, varcov_mCc,varXu,varYu, varXv, varYv)

    if ~opt.twodim
            varcov_mCw(c,:,:) = diag([varXu varYu]);
            
            if ~opt.white % if colored noise
                          den = varXu + varYu;
                        varXu = Puu(c)*varXu/den;
                        varYu = Puu(c)*varYu/den;
            varcov_mCc(c,:,:) = diag([varXu varYu]);
            end
            
                [sig1,sig2] = ut_linci(Xu(c),Yu(c),sqrt(varXu),sqrt(varYu));
            
            coef.Std.A(c,1) = sig1;  % Added (SI)
            coef.A_ci(c)    = 1.96*sig1;
            coef.Std.g(c,1) = sig2; % Added (SI)
            coef.g_ci(c)    = 1.96*sig2;
            
        else
            
            varcov_mCw(c,:,:) = diag([varXu varYu varXv varYv]);
            
            if ~opt.white
                  den = varXv + varYv;
                varXv = Pvv(c)*varXv/den;
                varYv = Pvv(c)*varYv/den;
                varcov_mCc(c,:,:) = diag([varXu varYu varXv varYv]);
            end
            
                 [sig1,sig2] = ut_linci(Xu(c)+1i*Xv(c),Yu(c)+1i*Yv(c),sqrt(varXu)+1i*sqrt(varXv),sqrt(varYu)+1i*sqrt(varYv));
            
            coef.Std.Lsmaj(c,1) = real(sig1); % Added (SI)
            coef.Std.Lsmin(c,1) = imag(sig1); % Added (SI)
            coef.Std.g(c,1)     = real(sig2); % Added (SI)
            coef.Std.theta(c,1) = imag(sig2); % Added (SI)
            
            coef.Lsmaj_ci(c)  = 1.96*real(sig1);
            coef.Lsmin_ci(c)  = 1.96*imag(sig1);
                coef.g_ci(c)  = 1.96*real(sig2);
            coef.theta_ci(c)  = 1.96*imag(sig2);
    end
end