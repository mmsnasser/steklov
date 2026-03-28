clear;clc;
%eig of e and k
n=2^10;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
et   = 8+5*exp(1i*t)+0.5*exp(6i*t);
etp  = 5i*exp(1i*t)+3i*exp(6i*t);
etpp = -5*exp(1i*t)-18*exp(6i*t);
Area = h*sum(real(et).*imag(etp));
alpha = 8;
%%
% E
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
% K
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
eige = eig(E); %eige = sort(eige);
eigk = eig(K); %eigk = sort(eigk);
%%
errk = NaN(size(eigk));
ind = abs(eigk)<1e-6; errk(ind)=abs(eigk(ind));
ind = abs(imag(eigk)-1)<1e-6; errk(ind)=abs(eigk(ind)-i);
ind = abs(imag(eigk)+1)<1e-6; errk(ind)=abs(eigk(ind)+i);
% 
erre = NaN(size(eige));
ind = abs(eige)<1e-6; erre(ind)=abs(eige(ind));
ind = abs(imag(eige)-1)<1e-6; erre(ind)=abs(eige(ind)-i);
ind = abs(imag(eige)+1)<1e-6; erre(ind)=abs(eige(ind)+i);
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(real(eigk),imag(eigk),'ob','MarkerFaceColor','b','MarkerSize',4);
axis square
axis([-1e-12 1e-12 -1.5 1.5])
xlabel('${\rm Re}[\hat\lambda]$')
ylabel('${\rm Im}[\hat\lambda]$')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
print -depsc eigK_ex3
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(real(eige),imag(eige),'ob','MarkerFaceColor','b','MarkerSize',4);
axis square
axis([-1e-12 1e-12 -1.5 1.5])
xlabel('${\rm Re}[\hat\lambda]$')
ylabel('${\rm Im}[\hat\lambda]$')
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
print -depsc eigE_ex3 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
semilogy(1:n,errk,'-b','LineWidth',1.5);
hold on; box on
semilogy(1:n,erre,'-r','LineWidth',1.5);
axis square
axis([0 1024 1e-16 1e-12])
xlabel('$k$')
ylabel('$|\hat\lambda_{k}-\hat\lambda_{k,n}|$')
legend('The matrix K','The matrix E','Location','northwest')
legend boxoff
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
print -depsc eigerr_ex3
%%