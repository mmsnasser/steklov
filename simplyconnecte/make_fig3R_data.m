clear;clc
%eigenvalues of KD using fft
% Cassini oval shape, Table 4 of Numerical studeies of the
% Steklov eigenvalue problem via convformal mappings
%%
% nv =[20:10:190,200:100:1000]';
nv=[20:20:400]';
for kk=1:length(nv)
kk
    n =nv(kk);

h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
w1   =  exp(1i*t);
w1p  =  1i*exp(1i*t);
w1pp =  -exp(1i*t);
et   = i*(0.4*w1.*sqrt(2./(1.16-0.84*w1.^2)));
etp  = i*(0.4*w1p.*sqrt(2./(1.16-0.84*w1.^2))+0.336*(w1.^2).*w1p.*...
          (sqrt(2)./(1.16-.84*w1.^2).^(3/2)));
etpp = i*(0.4*w1pp.*sqrt(2./(1.16-0.84*w1.^2))+1.008*w1.*(w1p.^2).*...
          (sqrt(2)./(1.16-.84*w1.^2).^(3/2))+0.336*(w1.^2).*w1pp.*...
          (sqrt(2)./(1.16-.84*w1.^2).^(3/2))+0.84672*(w1.^3).*(w1p.^2).*...
          (sqrt(2)./(1.16-.84*w1.^2).^(5/2)));
%
Area = h*sum(real(et).*imag(etp));
%
alpha = 0;
%%
A        =  et-alpha;
Ap       =  etp;
%
h        =   2*pi/n;
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
MD  = matD(n);
P   = diag(1./abs(etp));
%
KD = P*MD*E;
%%
KDI = KD+eye(n);
cnd(kk,1) = cond(KDI);
tic
[V,D] = eigs(KDI,12,'sm','Tolerance',1e-15,'MaxIterations',100,'Display',1);
timm(kk,1) = toc;
Dvn    = real(diag(D)-1);
%%
Lamex = [
          0.65940364662702
          2.31834186534589
          2.36353536012764
          2.68207118425903
          3.65243257022089
          4.04248974321298
          5.00265186479023
          5.07684221779547
          6.26494516290538
          6.34729870501637];
%
Lamn = Dvn(3:12);
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

%%
fig1=figure;
hold on; box on;
crv = et; crv(end+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',2);
set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',14)
set(ax,'LooseInset',get(ax,'TightInset'))
% print(fig1, 'fig_4air.pdf', '-dpdf', '-fillpage');
%%