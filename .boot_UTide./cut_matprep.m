function [yraw,Basis,coef,infer,nc,t] = cut_matprep(tin,uin,vin,lat,cnstit,varargin)
% S. Innocenti adapted from UT_SOLV1()
% 2019/11
% Set-up ut_tide regression for HA of a single record. See comments for UT_SOLV1().
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% 
%
%% input checking/conditioning, options parsing
[nt,t,u,v,tref,lor,~,opt,tgd,uvgd] = cut_slvinit(tin,uin,vin,cnstit,varargin);

%% constituent selection, initial storage
             opt.cnstit = cnstit;
[nNR,nR,nI,cnstit,coef] = ut_cnstitsel(tref,opt.rmin/(24*lor),opt.cnstit,opt.infer);
      coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat);
          coef.aux.opt  = opt;
          coef.aux.lat  = lat;
		  coef.aux.rpds = []; % initialize Residual Power Spectral Density (band averaged). Added (SI)

%% basis function matrix
[B,~,Q,beta,E,F,U,V] = cut_BEFUV(t,tref,cnstit,opt,lat,nt,nNR,nR,lor);

%% output
yraw = u;
if opt.twodim
   yraw = complex(u,v);
end

Basis.B = B;
Basis.E = E;
Basis.F = F;
Basis.U = U;
Basis.V = V;

infer.Q = Q;
infer.Beta = beta;

nc.nNR = nNR;
nc.nR  = nR;
nc.nI  = nI;
end