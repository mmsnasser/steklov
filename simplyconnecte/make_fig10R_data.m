clear;clc
%eigenvalues of KD using E
% Num ex Fig 4, table 1 of The exterior Steklov problem of Euclidean
% domains
%%
nv=[20:10:190,200:100:1000]';
% nv=[20:20:400]';
for kk=1:length(nv)
    kk
    n =nv(kk);
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
et    = 1.5*cos(t)+0.7*cos(2*t)-0.4+i*( -1.5*sin(t)-0.3*cos(t));
etp   =-1.5*sin(t)-1.4*sin(2*t)    +i*( -1.5*cos(t)+0.3*sin(t));
etpp  =-1.5*cos(t)-2.8*cos(2*t)    +i*(1.5*sin(t)+0.3*cos(t));
Area  = -h*sum(real(et).*imag(etp))
alpha = 0;
%%
for k=1:n
    for j=1:n
        if (k==j)
            N(k,j)  =  h*(1/pi)*imag(0.5*(etpp(k)/etp(k)));
            tM(k,j) =  h*(1/pi)*real(0.5*(etpp(k)/etp(k)));
        else
            N(k,j)  =  h*(1/pi)*imag((etp(j))/((et(j)-et(k))));
            tM(k,j) =  h*((1/(2*pi))*cot((t(k)-t(j))/2)...
                        +(1/pi)*real((etp(j))/((et(j)-et(k)))));
        end
    end
end
%
for k=1:n
    for j=1:n
        if (((j-k)/2)-floor((j-k)/2) == 0)
            K(k,j) = 0;
        else
            K(k,j) = (2/n)*cot(pi*(k-j)/n);
        end
    end
end
%
M        =  -K+tM;
E        =  -inv(eye(n)-N)*M;
%%
MD = matD(n);
P  = diag(1./abs(etp));
KD = P*MD*E;
%%
KDI = KD+eye(n);
cnd(kk,1) = cond(KDI);
tic
[V,D] = eigs(KDI,12,'sm','Tolerance',1e-15,'MaxIterations',100,'Display',1);
timm(kk,1) = toc;
%%
Dv    = real(diag(D)-1);
%%
Lamex = [
         0.544677700560795
         0.570816994124016
          1.12953414359678
          1.30930577399346
          1.74564067269481
          1.82146960379857
          2.29287781627365
          2.44997484632459
          2.90350756743372
          3.00356393166912];
%
Lamn = Dv(3:12);
for j=1:length(Lamex)
    rerrk(kk,j) =  abs((Lamn(j)-Lamex(j))/Lamex(j));
end
%%
end
%%
format short g
[nv rerrk]
[nv timm cnd]
format long g
%%
