
%%--------------------------------------------------------- 
function [astro,ader] = ut_astron(jd)
% UT_ASTRON()
% calculate astronomical constants
% input
%   jd = time [datenum UTC] (1 x nt)
% outputs
%   astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
%   ader = matrix of derivatives of astro [cycles/day] (6 x nt)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (copy of t_astron.m from t_tide, Pawlowicz et al 2002)

d=jd(:)'-datenum(1899,12,31,12,0,0);
D=d/10000;
args=[ones(size(jd));
      d;
      D.*D;
      D.^3];
sc= [ 270.434164,13.1763965268,-0.0000850, 0.000000039];
hc= [ 279.696678, 0.9856473354, 0.00002267,0.000000000];
pc= [ 334.329556, 0.1114040803,-0.0007739,-0.00000026];
npc=[-259.183275, 0.0529539222,-0.0001557,-0.000000050];
ppc=[ 281.220844, 0.0000470684, 0.0000339, 0.000000070];
astro=rem( [sc;hc;pc;npc;ppc]*args./360.0 ,1);
tau=rem(jd(:)',1)+astro(2,:)-astro(1,:);
astro=[tau;astro];
dargs=[zeros(size(jd));
       ones(size(jd));
       2.0e-4.*D;
       3.0e-4.*D.*D];
ader=[sc;hc;pc;npc;ppc]*dargs./360.0;
dtau=1.0+ader(2,:)-ader(1,:);
ader=[dtau;ader];