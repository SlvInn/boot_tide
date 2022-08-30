%%--------------------------------------------------------- 
function [Lsmaj,Lsmin,theta,g] = ut_cs2cep(XY)
% UT_CS2CEP()
% compute current ellipse parameters from cosine-sine coefficients
% inputs 
%   two-dim case: XY = [Xu Yu Xv Yv] 4-column matrix 
%   one-dim case: XY = [Xu Yu] 2-column matrix 
%                      OR 4-column w/ Xv=Yv=zeros(size(Xu))
%      where: Xu,Yu are cosine, sine coeffs of u, & same for v
% outputs
%   two-dim case:
%     Lsmaj, Lsmin = column vectors [units of XY] (size of Xu)
%            theta = column vector [deg. ccw rel. +x-axis, 0-180] (size of Xu)
%                g = column vector [degrees, 0-360] (size of Xu)
%   one-dim case: same, where Lsmaj = A, and Lsmin and theta are zeros
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

Xu = XY(:,1);
Yu = XY(:,2);
if size(XY,2)>2
    Xv = XY(:,3);
    Yv = XY(:,4);
else
    Xv = zeros(size(Xu));
    Yv = zeros(size(Yu));
end

   ap = ((Xu+Yv)+1i*(Xv-Yu))/2;
   am = ((Xu-Yv)+1i*(Xv+Yu))/2;
   Ap = abs(ap);
   Am = abs(am);
Lsmaj = Ap+Am;
Lsmin = Ap-Am;

% compute phase angle of the complex matrix
 epsp = angle(ap)*180/pi;
 epsm = angle(am)*180/pi;
theta = mod((epsp+epsm)/2,180);
    g = mod(-epsp+theta,360);
