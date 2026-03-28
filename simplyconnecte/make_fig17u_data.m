clear;clc
% eigenvalues of KD using E
%%
r  = [1 2 5 10]; r = r(:);
nv = [20:10:190,200:100:1000]';
%%
for m=1:length(r)
    for kk=1:length(nv)
    [m kk]
    n =nv(kk);
%
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
F=[]; MF=[]; K=[];N=[];tM=[];M=[];E=[];P=[];KD=[];KDI=[];V=[];D=[];
%%
MD = matD(n);
%%

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
    I     =  h*sum(sqrt(cos(t).^2+r(m)^2*sin(t).^2));
    a     =  2*pi/I;   
    et    =  a*( cos(t)-r(m)*1i*sin(t));
    etp   =  a*(-sin(t)-r(m)*1i*cos(t));
    etpp  =  a*(-cos(t)+r(m)*1i*sin(t));
    %
% 
    for k=1:n
        for j=1:n
            if (k==j)
              N(k,j)  =  h*(1/pi)*imag(0.5*(etpp(k)/etp(k)));
              tM(k,j) =  h*(1/pi)*real(0.5*(etpp(k)/etp(k)));
        else
              N(k,j)  =  h*(1/pi)*imag(etp(j)/(et(j)-et(k)));
              tM(k,j) =  h*((1/(2*pi))*cot((t(k)-t(j))/2)...
                        +(1/pi)*real(etp(j)/(et(j)-et(k))));
            end
        end
    end
    %
    M      =  -K+tM;
    E      =  -inv(eye(n)-N)*M;
    %
    P      =  diag(1./abs(etp));
    KD     =  P*MD*E;
    KDI    =  KD+eye(n);
    cnd(kk,m) = cond(KDI);
% 
    tic
    [V,D] = eigs(KDI,12,'sm','Tolerance',1e-15,'MaxIterations',100,'Display',1);
    timm(kk,m) = toc;
    Dv    = real(diag(D)-1);
%%
    end
end
%%
format short g
timm
cnd
format long g
%%
