function D = matD(n)
% 9/3/2025
% compute matrix L
%
F = dftmtx(n);
% 
R = diag([1:n/2-1]/n);
J = fliplr(eye(n/2-1));
W = zeros(n,n);
W(2:n/2,2:n/2)=R;
W(n/2+2:n,n/2+2:n)=-J*R*J;
W = -i*W;
D =  F*W*F';
% 
end