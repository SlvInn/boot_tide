function m_boot = cut_boot_reg(wB, yboot, w, opt) 
% S. Innocenti. 2019/09 - modified 2021/06
% Apply the least square estimation on the model defined by the B basis
% function and with the given opt 

switch opt.method

    case 'ols'
            iva = ~isnan(yboot);
            m_boot = wB(iva,:)\yboot(iva);

    case {'cauchy','andrews','bisquare','fair','huber', 'logistic','talwar','welsch','ridge', 'lasso', 'enet'}

        if any(strcmpi(opt.method,{'ridge', 'lasso', 'enet'}))
            warning([opt.method ' bootstrap not implemented yet >> using known weights'])
            assert (~isempty(w), [opt.method ' bootstrap need known weights'] )
        end

        if ~isempty(w)  
            % use the knwon weights in the regression  
                iva = ~isnan(yboot) & ~isnan(w);
            m_boot = wB(iva,:)\(w(iva).*yboot(iva));
        else
            % re-estimate weights
                iva = ~isnan(yboot);
            m_boot = robustfit(wB(iva,:),ctranspose(yboot(iva)),opt.method,opt.tunconst,'off');
        end

end   
end
        