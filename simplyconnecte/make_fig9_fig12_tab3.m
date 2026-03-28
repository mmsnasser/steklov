clear;clc
%eigenvalues of KD using E
% Num ex Fig 4, table 1 of The exterior Steklov problem of Euclidean
% domains both bounded and unbounded domains
%%
n =2^10;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
etb    = 1.5*cos(t)+0.7*cos(2*t)-0.4+i*( 1.5*sin(t)-0.3*cos(t));
etpb   =-1.5*sin(t)-1.4*sin(2*t)    +i*( 1.5*cos(t)+0.3*sin(t));
etppb  =-1.5*cos(t)-2.8*cos(2*t)    +i*(-1.5*sin(t)+0.3*cos(t));
etu    = 1.5*cos(t)+0.7*cos(2*t)-0.4+i*( -1.5*sin(t)-0.3*cos(t));
etpu   =-1.5*sin(t)-1.4*sin(2*t)    +i*( -1.5*cos(t)+0.3*sin(t));
etppu  =-1.5*cos(t)-2.8*cos(2*t)    +i*(1.5*sin(t)+0.3*cos(t));
Area   = h*sum(real(etb).*imag(etpb));
alpha  = 0;
%%
fig1=figure;
hold on; box on;
crv = etb; crv(end+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',2);
plot(real(alpha),imag(alpha),'sr','MarkerFaceColor','r');
set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
print -depsc fig_ex5_e
%%
A        =  etb-alpha;
Ap       =  etpb;
% 
for k=1:n
    for j=1:n
        if (k==j)
            Nb(k,j)  =  h*(1/pi)*imag(0.5*(etppb(k)/etpb(k))-Ap(k)/A(k));
            tMb(k,j) =  h*(1/pi)*real(0.5*(etppb(k)/etpb(k))-Ap(k)/A(k));
        else
            Nb(k,j)  =  h*(1/pi)*imag((A(k)*etpb(j))/(A(j)*(etb(j)-etb(k))));
            tMb(k,j) =  h*((1/(2*pi))*cot((t(k)-t(j))/2)...
                        +(1/pi)*real((A(k)*etpb(j))/(A(j)*(etb(j)-etb(k)))));
        end
    end
end
%
for k=1:n
    for j=1:n
        if (k==j)
            Nu(k,j)  =  h*(1/pi)*imag(0.5*(etppu(k)/etpu(k)));
            tMu(k,j) =  h*(1/pi)*real(0.5*(etppu(k)/etpu(k)));
        else
            Nu(k,j)  =  h*(1/pi)*imag((etpu(j))/((etu(j)-etu(k))));
            tMu(k,j) =  h*((1/(2*pi))*cot((t(k)-t(j))/2)...
                        +(1/pi)*real((etpu(j))/((etu(j)-etu(k)))));
        end
    end
end
%
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
Mb        =  -K+tMb;
Eb        =  -inv(eye(n)-Nb)*Mb;
%
Mu        =  -K+tMu;
Eu        =  -inv(eye(n)-Nu)*Mu;
%%
MD = matD(n);
Pb  = diag(1./abs(etpb));
Pu  = diag(1./abs(etpu));
% 
KDb = Pb*MD*Eb;
KDIb = KDb+eye(n);
[Vb,Db] = eigs(KDIb,102,'smallestabs');
% 
KDu = Pu*MD*Eu;
KDIu = KDu+eye(n);
[Vu,Du] = eigs(KDIu,102,'smallestabs');
%%
format long g
Dvb    = real(diag(Db)-1);
Dvu    = real(diag(Du)-1);
[Dvb(3:12) Dvu(3:12)]
%%
L   =  h*sum(abs(etpb));
%%
v= 1:100; v= v';
xeb = Dvb(3:102);
xeu = Dvu(3:102);
%
% for k=1:50
%    lam(2*k-1) = (2*pi*k)/L;
%    lam(2*k)   = (2*pi*k)/L;
% end
% lam=lam';
%
figure;
hold on; box on;
plot(v,xeb,'b','LineWidth',2);
plot(v,xeu,'r','LineWidth',2);
plot(v,(pi*v)/L,':k','LineWidth',2);
% plot(v,lam,':k','LineWidth',2);
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
legend({'$\lambda_k$ for $G_1$','$\lambda_k$ for $G_2$','$\frac{\pi k}{|\Gamma|}$'},...
    'Location','northwest','Interpreter','latex')
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',14)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
%%
figure;
hold on; box on;
plot(v(1:10),xeb(1:10),'^b','LineWidth',2,'MarkerFaceColor','b');
plot(v(1:10),xeu(1:10),'vr','LineWidth',2,'MarkerFaceColor','r');
plot(v(1:10),(pi*v(1:10))/L,':k','LineWidth',2);
% plot(v,lam,':k','LineWidth',2);
set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
legend({'$G_1$ ','$G_2$ ','$\frac{\pi k}{|\Gamma|}$'},...
    'Location','northwest','Interpreter','latex')
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_ex5_10eig_1
%%
figure;
hold on; box on;
plot(v(30:40),xeb(30:40),'^b','LineWidth',2,'MarkerFaceColor','b');
plot(v(30:40),xeu(30:40),'vr','LineWidth',2,'MarkerFaceColor','r');
plot(v(30:40),(pi*v(30:40))/L,':k','LineWidth',2);
% plot(v,lam,':k','LineWidth',2);
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
legend({'$G_1$ ','$G_2$ ','$\frac{\pi k}{|\Gamma|}$'},...
    'Location','northwest','Interpreter','latex')
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
xticks([30:2:40])
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_ex5_10eig_2
%%
figure;
hold on; box on;
plot(v(90:100),xeb(90:100),'^b','LineWidth',2,'MarkerFaceColor','b');
plot(v(90:100),xeu(90:100),'vr','LineWidth',2,'MarkerFaceColor','r');
plot(v(90:100),(pi*v(90:100))/L,':k','LineWidth',2);
% plot(v,lam,':k','LineWidth',2);
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
legend({'$G_1$ ','$G_2$ ','$\frac{\pi k}{|\Gamma|}$'},...
    'Location','northwest','Interpreter','latex')
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_ex5_10eig_3
%%