
%%--------------------------------------------------------- 
function ain = ut_cluster(ain,clusang)
% UT_CLUSTER()
% UTide v1p0 9/2011 d.codiga@gso.uri.edu
% (copy of t_cluster.m of t_tide, Pawlowicz et al 2002)
ii=(ain-ain(:,ones(1,size(ain,2))))>clusang/2;
ain(ii)=ain(ii)-clusang;
ii=(ain-ain(:,ones(1,size(ain,2))))<-clusang/2;
ain(ii)=ain(ii)+clusang;
end
