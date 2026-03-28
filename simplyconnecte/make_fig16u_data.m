clear;clc
% eigenvalues of KD using E
%%
n =2^11;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
MD = matD(n);
%%
r = [1;2;5;10];
rk = length(r);
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
alpha = 0;
%
for m=1:rk
    m
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
    [V,D]  =  eigs(KDI,12,'smallestabs');
    Dv    = real(diag(D)-1);
%%
for j=1:10
        xe(j,m) = Dv(j+2);
end
%%
end
%%
