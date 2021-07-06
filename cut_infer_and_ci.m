function [coef] = cut_infer_and_ci(coef,varcov_mCw,varcov_mCc,nNR,nI)
% S. Innocenti adapted from ut_solv1() [UTide v1p0 9/2011 d.codiga@gso.uri.edu]
% 2019/09

% complex coeffs
    apI = nan*ones(1,nI);
    amI = apI;
    ind = 1;
    
for k = 1:nR
    apI(ind:ind+cnstit.R{k}.nI-1) = cnstit.R{k}.I.Rp*ap(nNR+k);
    amI(ind:ind+cnstit.R{k}.nI-1) = cnstit.R{k}.I.Rm*am(nNR+k);
                              ind = ind + cnstit.R{k}.nI;
end   

% cos/sin coeffs and c.e.p.
XuI = real(apI+amI)';
YuI = -imag(apI-amI)';
    
if ~opt.twodim
    [A,~,~,g] = ut_cs2cep([XuI YuI]);
        coef.A = [coef.A; A;];
        coef.g = [coef.g; g;];
else
                        XvI = imag(apI+amI)';
                        YvI = real(apI-amI)';
    [Lsmaj,Lsmin,theta,g] = ut_cs2cep([XuI YuI XvI YvI]);
                coef.Lsmaj = [coef.Lsmaj; Lsmaj;];
                coef.Lsmin = [coef.Lsmin; Lsmin;];
                coef.theta = [coef.theta; theta;];
                    coef.g = [coef.g; g;];
end
    
% confidence intervals
if opt.linci % linearized
    for k = 1:nR
        if opt.white
            if ~opt.twodim
                varReap = 0.25*varcov_mCw(nNR+k,1,1);
                varImap = 0.25*varcov_mCw(nNR+k,2,2);
            else
                varReap = 0.25*(varcov_mCw(nNR+k,1,1) + varcov_mCw(nNR+k,4,4));
                varImap = 0.25*(varcov_mCw(nNR+k,2,2) + varcov_mCw(nNR+k,3,3));
            end
        else
            if ~opt.twodim
                varReap = 0.25*varcov_mCc(nNR+k,1,1);
                varImap = 0.25*varcov_mCc(nNR+k,2,2);
            else
                varReap = 0.25*(varcov_mCc(nNR+k,1,1) + varcov_mCc(nNR+k,4,4));
                varImap = 0.25*(varcov_mCc(nNR+k,2,2) + varcov_mCc(nNR+k,3,3));
            end
        end
        varXuHH = (real(cnstit.R{k}.I.Rp).^2 + real(cnstit.R{k}.I.Rm).^2)*varReap + ...
                    (imag(cnstit.R{k}.I.Rp).^2 + imag(cnstit.R{k}.I.Rm).^2)*varImap;
        varYuHH = (imag(cnstit.R{k}.I.Rp).^2 + imag(cnstit.R{k}.I.Rm).^2)*varReap + ...
                    (real(cnstit.R{k}.I.Rp).^2 + real(cnstit.R{k}.I.Rm).^2)*varImap;
        for c = 1 : length(varXuHH)
            if ~opt.twodim
                [sig1,sig2] = ut_linci(Xu(nNR+k),Yu(nNR+k),sqrt(varXuHH(c)),sqrt(varYuHH(c)));
                    coef.A_ci = [coef.A_ci; 1.96*sig1;];
                    coef.g_ci = [coef.g_ci; 1.96*sig2;];
            else
                [sig1,sig2]= ut_linci(Xu(nNR+k)+1i*Xv(nNR+k),Yu(nNR+k)+1i*Yv(nNR+k),...
                                sqrt(varXuHH(c))+1i*sqrt(varYuHH(c)),...
                                sqrt(varYuHH(c))+1i*sqrt(varXuHH(c)));
                coef.Lsmaj_ci = [coef.Lsmaj_ci; 1.96*real(sig1);];
                coef.Lsmin_ci = [coef.Lsmin_ci; 1.96*imag(sig1);];
                coef.g_ci     = [coef.g_ci; 1.96*real(sig2);];
                coef.theta_ci = [coef.theta_ci; 1.96*imag(sig2);];
            end
        end
    end
    

else % monte carlo
    ind = 0;
    for k = 1:nR
        if ~opt.twodim
            if opt.white
                mCall = mvnrnd([Xu(nNR+k) Yu(nNR+k)], squeeze(varcov_mCw(nNR+k,:,:)),opt.nrlzn);
            else
                mCall = mvnrnd([Xu(nNR+k) Yu(nNR+k)],squeeze(varcov_mCc(nNR+k,:,:)),opt.nrlzn);
            end
            
            
            [A,~,~,g] = ut_cs2cep(mCall);
                    g(1) = coef.g(nNR+k);
                    g = ut_cluster(g,360);
                    
            for lk = 1 : cnstit.R{k}.nI
                ind = ind+1;
                        apI = 0.5*cnstit.R{k}.I.Rp(lk)*(A.*exp(-1i*g*pi/180));
                        amI = conj(apI);
                        XuI = real(apI+amI);
                        YuI = -imag(apI-amI);
                    [A,~,~,g] = ut_cs2cep([XuI YuI]);
                    
                        A_ci = median(abs(A-median(A)))/0.6745;
                coef.A_ci = [coef.A_ci; 1.96*A_ci;];
                
                coef.A_rlzn = [coef.A_rlzn; A';]; %(PM): added 
                            g = ut_cluster(g,360);
                            
                        g_ci = median(abs(g-median(g)))/0.6745;
                coef.g_rlzn = [coef.g_rlzn; g';]; %(PM): added 
                    coef.g_ci = [coef.g_ci; 1.96*g_ci;];
            end
        else
            if opt.white
                mCall = mvnrnd([Xu(nNR+k) Yu(nNR+k) Xv(nNR+k) Yv(nNR+k)],squeeze(varcov_mCw(nNR+k,:,:)), opt.nrlzn);
            else
                mCall = mvnrnd([Xu(nNR+k) Yu(nNR+k) Xv(nNR+k) Yv(nNR+k)],squeeze(varcov_mCc(nNR+k,:,:)), opt.nrlzn);
            end                    
            [Lsmaj,Lsmin,theta,g] = ut_cs2cep(mCall);
                             g(1) = coef.g(nNR+k);
                                g = ut_cluster(g,360);
                         theta(1) = coef.theta(nNR+k);
                            theta = ut_cluster(theta,360);
                            
            for lk = 1:cnstit.R{k}.nI
                
                ind = ind+1;
                apI = cnstit.R{k}.I.Rp(lk).*0.5*(Lsmaj+Lsmin).* exp(1i*(theta-g)*pi/180);
                amI = cnstit.R{k}.I.Rm(lk).*0.5*(Lsmaj-Lsmin).* exp(1i*(theta+g)*pi/180);
                XuI = real(apI+amI);
                YuI = -imag(apI-amI);
                XvI = imag(apI+amI);
                YvI = real(apI-amI);
                
                [Lsmaj,Lsmin,theta,g] = ut_cs2cep([XuI YuI XvI YvI]);
                             Lsmaj_ci = median(abs(Lsmaj-median(Lsmaj)))/0.6745;
                        coef.Lsmaj_ci = [coef.Lsmaj_ci; 1.96*Lsmaj_ci;];
                             Lsmin_ci = median(abs(Lsmin-median(Lsmin)))/0.6745;
                        coef.Lsmin_ci = [coef.Lsmin_ci; 1.96*Lsmin_ci;];
                                theta = ut_cluster(theta,360);
                             theta_ci = median(abs(theta-median(theta)))/0.6745;
                        coef.theta_ci = [coef.theta_ci; 1.96*theta_ci;];
                                    g = ut_cluster(g,360);
                                 g_ci = median(abs(g-median(g)))/0.6745;
                            coef.g_ci = [coef.g_ci; 1.96*g_ci;];
            end
        end            
    end
end
    
