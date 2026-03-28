clear;clc
%eigenvalues of KD using fft
% Table 3 of Numerical studeies of the
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
et   = 8+5*exp(i*t)+0.5*exp(6i*t);
etp  = 5i*exp(1i*t)+3i*exp(6i*t);
etpp = -5*exp(1i*t)-18*exp(6i*t);
Area = h*sum(real(et).*imag(etp));
%
alpha = 8;
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
MD = matD(n);
P  = diag(1./abs(etp));
KD = P*MD*E;
%%
KDI = KD+eye(n);
cnd(kk,1) = cond(KDI);
tic
[V,D] = eigs(KDI,12,'sm','Tolerance',1e-15,'MaxIterations',100,'Display',1);
timm(kk,1) = toc;
Dvn    = real(diag(D)-1);
%%
format long g
Dvn(1:12)
%%
Lamex = [
         0.176962409070052
         0.176962409070062
         0.326314227252894
          0.32631422725292
         0.600966675036201
         0.600966675036212
          0.73515390643172
         0.735153906431743
         0.839271996894644
         0.988527105824377];
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