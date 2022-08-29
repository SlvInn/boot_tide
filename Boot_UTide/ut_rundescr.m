
%%--------------------------------------------------------- 
function rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat)
% UT_RUNDESCR()
% build character array with descriptions of run
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

       rd = cell(10+nI,1);
    rd{1} = sprintf('UTide Run Description (executed %s):',datestr(now));
    if opt.equi
        if sum(tgd)>sum(uvgd)
            rd{2} = sprintf(['  Input(%dD): %d values (gappy); %d '...
                'equisp. times, dt=%#.2g hr.'],...
                opt.twodim+1,sum(uvgd),sum(tgd),24*mode(diff(t)));
        else
            rd{2} = sprintf(['  Input(%dD): %d values, equispaced times,'...
                ' dt=%#.3g hr.'],...
                opt.twodim+1,sum(uvgd),24*mode(diff(t)));
        end
    else
        rd{2} = sprintf(['  Input(%dD): %d vals; irreg., dt(min/mean/max)='...
            ' %#.2g/%#.2g/%#.2g hr.'],opt.twodim+1,length(t),...
            24*min(diff(t)),24*mean(diff(t)),24*max(diff(t)));
    end
    if max(t)-min(t)<365.25
        rd{3} = sprintf('    dur''n %#.3g dy; %s to %s GMT.',...
            max(t)-min(t),datestr(min(t),'dd-mmm-yy HH:MM'),...
            datestr(max(t),'dd-mmm-yy HH:MM'));
    else
        rd{3} = sprintf('    dur''n %#.3g yr; %s to %s GMT.',...
            (max(t)-min(t))/365.25,datestr(min(t),'dd-mmm-yy HH:MM'),...
            datestr(max(t),'dd-mmm-yy HH:MM'));
    end
    if isequal(opt.method,'ols')
        rd{4} = '  Method: OLS.';
    else
        rd{4} = sprintf(['  Method: IRLS, %s wt func; tunprm red''n fac '...
            '%#.3g.'],opt.method,opt.tunrdn);
    end
    if opt.linci
        rd{5} = '  Conf-ints: linearized;';
    else
        rd{5} = sprintf('  Conf-ints: Monte-Carlo, %d realizations;',...
            opt.nrlzn);
    end
    if opt.white
        rd{5} = [rd{5} ' white.'];
    else
        rd{5} = [rd{5} ' colored'];
        if ~opt.equi
            rd{5} = [rd{5} sprintf(' (Lmb-Scg ovrsmp %d).',...
                opt.lsfrqosmp)];
        else
            rd{5} = [rd{5} '.'];
        end
    end
    rd{6} = sprintf('  Model: allc=%d cnstits (%d non-ref), ',nNR+nR+nI,nNR);
    if isequal(opt.cnstit,'auto')
        rd{6} = [rd{6} sprintf('auto (Rmin=%#.3g);',...
            opt.rmin)];
    else
        rd{6} = [rd{6} 'specified as inputs;'];
    end
    if opt.notrend
        rd{7} = '    no trend included;';
    else
        rd{7} = '    trend included;';
    end
    if isempty(opt.prefilt)
        rd{7} = [rd{7} ' no prefiltering correction;'];
    else
        rd{7} = [rd{7} ' prefiltering correction applied;'];
    end
    if opt.nodsatnone
        rd{8} = '    no nod/sat corrections;';
    elseif opt.nodsatlint
        rd{8} = sprintf(['    nod/sat corrections (lat=%.3f), '...
            'linearized times;'],lat);
    else
        rd{8} = sprintf(['    exact nod/sat corrections (lat=%.3f);'],lat);
    end
    if opt.gwchnone
        rd{9} = '    g= raw phases (no Gwich correction);';
    elseif opt.gwchlint
        rd{9} = '    g= Gwich phase lag, astr arg linearized times;';
    else
        rd{9} = '    g= Gwich phase lag, astr arg exact times;';
    end
    if isempty(opt.infer)
        rd{10} = '    no inference.';
    else
        rd{10} = '    exact inference:';
        if opt.inferaprx
            rd{10} = '    approximate inference:';
        end
        for j = 1:nI
            if opt.twodim
                rd{10+j} = sprintf(['     %s inferred from %s r+,r-= '...
                    '%#.3g,%#.3g zta+,zta-= %#.3g,%#.3g deg.'],...
                    deblank(opt.infer.infnam{j}),...
                    deblank(opt.infer.refnam{j}),...
                    opt.infer.amprat(j),opt.infer.amprat(nI+j),...
                    opt.infer.phsoff(j),opt.infer.phsoff(nI+j));
            else
                rd{10+j} = sprintf(['     %s inferred from %s r= '...
                    '%#.3g zta= %#.3g deg.'],...
                    deblank(opt.infer.infnam{j}),...
                    deblank(opt.infer.refnam{j}),...
                    opt.infer.amprat(j),opt.infer.phsoff(j));
            end
        end
    end
    rundescr = char(rd);
end
