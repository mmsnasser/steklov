clear;clc
% 
%%
riv = [1.0;1.9838737083599;3.1178117418790;4.2783368585899;5.44203249398502;...
           6.6045610453282;7.7653565042885;8.9245669181406;10.0824459699803];
%
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10,'MaxIter',100);
% 
resm =[];
for k=1:length(riv)
    % 
    f = @(r) lambdk(k,r);
    ri = riv(k);
    r= fminsearch(f,ri,options)
    %
    Lamr2 = lambdk(k,r,10)
    %
    resv = []; resv = [r; Lamr2];
    resm = [resm resv];
end
%
%%
resm
%%