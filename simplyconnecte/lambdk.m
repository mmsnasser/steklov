function dL = lambdk(Lk,r,nk)
% find the value of lambda_(k+1)-lambda_(k) for a given value of r
%%
if r<=5
    n =2^10;
else
    n=2^11;
end
%
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%
MD = matD(n);
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
alpha = 0;
%
I     =  h*sum(sqrt(cos(t).^2+r^2*sin(t).^2));
a     =  2*pi/I;   
et    =  a*(cos(t)+r*1i*sin(t));
etp   =  a*(-sin(t)+r*1i*cos(t));
etpp  =  a*(-cos(t)-r*1i*sin(t));
%
A        =  et-alpha;
Ap       =  etp;
%  
for k=1:n
    for j=1:n
        if (k==j)
            N(k,j)  =  h*(1/pi)*imag(0.5*(etpp(k)/etp(k))-Ap(k)/A(k));
            tM(k,j) =  h*(1/pi)*real(0.5*(etpp(k)/etp(k))-Ap(k)/A(k));
        else
            N(k,j)  =  h*(1/pi)*imag((A(k)*etp(j))/(A(j)*(et(j)-et(k))));
            tM(k,j) =  h*((1/(2*pi))*cot((t(k)-t(j))/2)...
                +(1/pi)*real((A(k)*etp(j))/(A(j)*(et(j)-et(k)))));
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
[V,D]  =  eigs(KDI,102,'smallestabs');
Dv     =  real(diag(D)-1);
%
if nargin==2
    dL = Dv(Lk+3)-Dv(Lk+2);
elseif nargin==3
    dL = Dv(3:nk+2);
end
%% 
end