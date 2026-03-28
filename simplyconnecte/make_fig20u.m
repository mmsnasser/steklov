clear;clc
% eigenvalues of KD using E
%%
r = 1:5; r = r(:);
rk = length(r);
%%
for m=1:rk
    m
    n =2^10;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
MD = matD(n);
%%
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
    I     =  h*sum(sqrt(cos(t).^2+r(m)^2*sin(t).^2));
    a     =  2*pi/I;   
    et    =  a*(cos(t)-r(m)*1i*sin(t));
    etp   =  a*(-sin(t)-r(m)*1i*cos(t));
    etpp  =  a*(-cos(t)+r(m)*1i*sin(t));
    %
    Area = -h*sum(real(et).*imag(etp));
    xa(m) = 1/(a.*sqrt(r(m)));
    %
%%
    
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
    M        =  -K+tM;
    E        =  -inv(eye(n)-N)*M;
    P  = diag(1./abs(etp));
    %
    KD     =  P*MD*E;
    KDI    =  KD+eye(n);
    [V,D]  =  eigs(KDI,102,'smallestabs');
    Dv    = real(diag(D)-1);
%%
for j=1:100
        xe(j,m) = Dv(j+2);
end
%%
end
%%
clr = ['m','c','g','b','r'];
mrk = ['d','s','^','v','p'];
v= 1:100; v= v';
figure;
hold on; box on
for mm=1:length(r)
plot(v,xe(:,r(mm)),'color',clr(mm),'LineWidth',1.5);
end
plot(v,v/2,':k','LineWidth',1.5)
legend({'$r=1$','$r=2$','$r=3$','$r=4$','$r=5$','$k/2$'},...
    'Location','northwest','Interpreter','latex','NumColumns',2)
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
axis([0 100 0 55])
xticks([0:20:100])
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
%%
%%
figure(11);clf
hold on; box on
for mm=1:length(r)
plot(v(1:10),xe(1:10,r(mm)),'Marker',mrk(mm),'LineWidth',1.5,'LineStyle','none');
end
plot(v(1:10),v(1:10)/2,':k','LineWidth',1.5)
legend({'$r=1$','$r=2$','$r=3$','$r=4$','$r=5$','$k/2$'},...
    'Location','northwest','Interpreter','latex','NumColumns',2)
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
axis([0 10 0 6])
xticks([0:2:10])
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_ex_ell_eig10_2_U
%%
figure(12);clf
hold on; box on
for mm=1:length(r)
plot(v(90:100),xe(90:100,r(mm)),'Marker',mrk(mm),'LineWidth',1.5,'LineStyle','none');
end
plot(v(90:100),v(90:100)/2,':k','LineWidth',1.5)
legend({'$r=1$','$r=2$','$r=3$','$r=4$','$r=5$','$k/2$'},...
    'Location','northwest','Interpreter','latex','NumColumns',2)
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
axis([90 100 45 51])
xticks([90:2:100])
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_ex_ell_eig10_3_U
%%
figure
hold on; box on
for mm=1:length(r)
plot(v(30:40),xe(30:40,r(mm)),'Marker',mrk(mm),'LineWidth',1.5,'LineStyle','none');
end
plot(v(30:40),v(30:40)/2,':k','LineWidth',1.5)
legend({'$r=1$','$r=2$','$r=3$','$r=4$','$r=5$','$k/2$'},...
    'Location','northwest','Interpreter','latex','NumColumns',2)
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
% axis([90 100 45 51])
xticks([30:2:40])
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$k$','Interpreter','latex')
ylabel('$\lambda_k$','Interpreter','latex')
print -depsc fig_ex_ell_eig10_4_U
%%