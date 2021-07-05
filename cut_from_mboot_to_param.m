function coef_boot = cut_from_mboot_to_param(coef_boot,n_boot,opt,nNR,nR) %add Q,beta,cnstit,lor for inference
% Pass from complex or sin-cos regression parameters to 
% tidal amplitude and phases for each bootstrap resample

    for bs = 1 : n_boot
        m_boot = coef_boot.M(:,bs);

        if opt.complexf==1

            % from Complex to Amp/Pha coeff 
            ap = m_boot([1:nNR       2*nNR+(1:nR)]);
            am = m_boot([nNR+(1:nNR) 2*nNR+nR+(1:nR)]);
            if ~isempty(opt.infer)
                if opt.inferaprx
                    for k = 1:nR
                            ap(nNR+k) = ap(nNR+k)./(1+beta(k)*cnstit.R{k}.I.Rp.*Q(k));
                            am(nNR+k) = am(nNR+k)./(1+beta(k)*cnstit.R{k}.I.Rm.*conj(Q(k)));
                    end
                end
            end
            Xu = real(ap+am);
            Yu = -imag(ap-am);
            if ~opt.twodim
                [coef_boot.A(:,bs),~,~,coef_boot.g(:,bs)] = ut_cs2cep([Xu Yu]);
            else
                Xv = imag(ap+am);
                Yv = real(ap-am);

                [coef_boot.Lsmaj(:,bs),coef_boot.Lsmin(:,bs),coef_boot.theta(:,bs),coef_boot.g(:,bs)] = ut_cs2cep([Xu Yu Xv Yv]);
            end

   
            % mean and trend
            if opt.twodim
                if opt.notrend
                    coef_boot.umean(1,bs) = real(m_boot(end));
                    coef_boot.vmean(1,bs) = imag(m_boot(end));
                else
                    coef_boot.umean(1,bs) = real(m_boot(end-1));
                    coef_boot.vmean(1,bs) = imag(m_boot(end-1));
                    coef_boot.uslope(1,bs) = real(m_boot(end))/lor;
                    coef_boot.vslope(1,bs) = imag(m_boot(end))/lor;
                end
            else
                if opt.notrend
                    coef_boot.mean(1,bs) = real(m_boot(end));
                else
                    coef_boot.mean(1,bs) = real(m_boot(end-1));
                    coef_boot.slope(1,bs) = real(m_boot(end))/lor;
                end
            end

        else
            
            
            ap = m_boot([1:nNR       2*nNR+(1:nR)]);
            am = m_boot([nNR+(1:nNR) 2*nNR+nR+(1:nR)]);
            if ~isempty(opt.infer)
                error 'code for the opt.infer option must be done: see cut_fromComplex2AmpPha.m ~lines 6-13'
                % Q,beta, and cnstit would be used here
            end

            coef_boot.A(:,bs) = sqrt(ap.^2 + am.^2);
                          pha = atan2d(am,ap);
                          pha = ut_cluster(pha,360);
                   pha(pha<0) = 360+pha(pha<0);
            coef_boot.g(:,bs) = pha;



            %% mean and trend
            if opt.notrend
                coef_boot.mean(1,bs)  = real(m_boot(end));
            else
                coef_boot.mean(1,bs)  = real(m_boot(end-1));
                coef_boot.slope(1,bs) = real(m_boot(end))/lor;
            end


        end

    end
end    

