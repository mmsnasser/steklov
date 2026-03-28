clear;clc
%eigenvalues of KD using fft
% Cassini oval shape, Table 4 of Numerical studeies of the
% Steklov eigenvalue problem via convformal mappings
%%
n =2^10;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
% syms f(w)
% f   = 0.4*w.*sqrt(2/(1.16-0.84*w^2));
% fp  = diff(f,1);
% fpp = diff(f,2);
%
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
% et   = double(subs(f,exp(i*t)));
% etp  = double(subs(fp,exp(i*t)))*i.*exp(i*t);
% etpp = double(subs(fpp,exp(i*t)).*exp(i*t)+subs(fp,exp(i*t)))*(-1).*exp(i*t);
% 
Area = h*sum(real(et).*imag(etp));
L    = h*sum(abs(etp));
alpha = 0;
%%
fig1=figure;
hold on; box on;
crv = et; crv(end+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',2);
plot(real(alpha),imag(alpha),'sr','MarkerFaceColor','r');
set(fig1,'PaperSize',[5  5]);
grid on; axis square
axis([-1 1 -1 1])
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
% print(fig1, 'fig_4air.pdf', '-dpdf', '-fillpage');
print -depsc fig_ex4e
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
format long g
%
KDI = KD+eye(n);
[V,D] = eigs(KDI,102,'smallestabs');
Dv    = real(diag(D)-1)*sqrt(Area);
lam_ns = real(diag(D)-1);
Dv(1:12)
lam_ns(1:12)
%%
v= 1:100; v= v';
xe = lam_ns(3:102);
%
figure;
hold on; box on;
plot(v,xe,'b','LineWidth',2);
plot(v,(pi*v)/L,':k','LineWidth',2);
set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
% print -depsc fig_ex4e_100eig
%%
figure;
hold on; box on;
plot(v(1:10),xe(1:10),'ob','LineWidth',2,'MarkerFaceColor','b');
plot(v(1:10),(pi*v(1:10))/L,':k','LineWidth',2);
% plot(v,lam,':k','LineWidth',2);
grid on; axis square
ax=gca; 
legend({'$\lambda_k$','$\frac{\pi k}{|\Gamma|}$'},...
    'Location','northwest','Interpreter','latex')
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_4_10eig_1
%%
figure;
hold on; box on;
plot(v(30:40),xe(30:40),'ob','LineWidth',2,'MarkerFaceColor','b');
plot(v(30:40),(pi*v(30:40))/L,':k','LineWidth',2);
% plot(v,lam,':k','LineWidth',2);
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
legend({'$\lambda_k$','$\frac{\pi k}{|\Gamma|}$'},...
    'Location','northwest','Interpreter','latex')
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
xticks([30:2:40])
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_4_10eig_2
%%
figure;
hold on; box on;
plot(v(90:100),xe(90:100),'ob','LineWidth',2,'MarkerFaceColor','b');
plot(v(90:100),(pi*v(90:100))/L,':k','LineWidth',2);
% plot(v,lam,':k','LineWidth',2);
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
legend({'$\lambda_k$','$\frac{\pi k}{|\Gamma|}$'},...
    'Location','northwest','Interpreter','latex')
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_4_10eig_3
%%