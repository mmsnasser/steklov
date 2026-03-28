clear;clc
% eigenvalues of KD using E
%%
rv = [1 2 5 10]';
nv = [20:20:400]';
%
xexact = [
         0.999999999999933          0.58949168152398         0.270355907931518         0.140617215789966
         0.999999999999947          1.67277042661387         0.905110688529072         0.489511425425261
          1.99999999999996          1.68973486925529          1.80635304751159           1.0280361683962
          1.99999999999997          2.38532840049433          2.86947989369561          1.73305572745026
          2.99999999999996          2.81816245473333          3.57053142746906          2.57726697480456
          2.99999999999998          3.20265984655634           4.0071964304634          3.53181196769968
          3.99999999999995          3.90690475610796          4.16549081354628          4.56871900960875
          3.99999999999997          4.10374347987443           4.8180057759941          5.66278638318708
          4.99999999999992          4.95199807234296          5.16214524419659          6.71145422035504
          4.99999999999997          5.05312773715011          5.52368346740819          6.79275482904478
          ];
%
%%
for m=1:length(rv)
    r = rv(m);
    for kk=1:length(nv)
        n = nv(kk);
        [r n]
%
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
F=[]; MF=[]; K=[];N=[];tM=[];M=[];E=[];P=[];KD=[];KDI=[];V=[];D=[];
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

    I     =  h*sum(sqrt(cos(t).^2+r^2*sin(t).^2));
    a     =  2*pi/I;   
    et    =  a*( cos(t)+r*1i*sin(t));
    etp   =  a*(-sin(t)+r*1i*cos(t));
    etpp  =  a*(-cos(t)-r*1i*sin(t));
    %
    A        =  et-alpha;
    Ap       =  etp;
    %
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
    MD = matD(n);
    P      =  diag(1./abs(etp));
    KD     =  P*MD*E;
    KDI    =  KD+eye(n);
% 
    [V,D]  =  eigs(KDI,12,'smallestabs');
    Dv    = real(diag(D)-1);
    %
    for j=1:10
        xe(j,m) = Dv(j+2);
    end
    %
    for j=1:10
        reerr{m}(kk,j) =  abs((xexact(j,m)-xe(j,m))/xexact(j,m));
    end
    end
    %
    reerr{m}(reerr{m}==0)=eps;
end
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
markers = {'o', '+', '*', 'x', 'square','diamond','^','v','<','>','pentagram','hexagram','-','|','.'};
%
for j=1:10
    semilogy(nv,reerr{1}(:,j),'Marker',markers(j),'LineWidth',1.5)
    hold on; box on
end
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Relative Error','FontSize',14,'Interpreter','latex');
legend({'$\lambda_1$','$\lambda_2$','$\lambda_3$','$\lambda_4$',...
    '$\lambda_5$','$\lambda_6$','$\lambda_7$','$\lambda_8$','$\lambda_9$',...
    '$\lambda_{10}$'},...
    'Location','northeast','Interpreter','latex','NumColumns',3)
axis([0 400 1e-15  1e0])
xticks([0:100:400])
yticks([10.^[-15:5:0]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ell_err10_1
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
markers = {'o', '+', '*', 'x', 'square','diamond','^','v','<','>','pentagram','hexagram','-','|','.'};
%
for j=1:10
    semilogy(nv,reerr{2}(:,j),'Marker',markers(j),'LineWidth',1.5)
    hold on; box on
end
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Relative Error','FontSize',14,'Interpreter','latex');
legend({'$\lambda_1$','$\lambda_2$','$\lambda_3$','$\lambda_4$',...
    '$\lambda_5$','$\lambda_6$','$\lambda_7$','$\lambda_8$','$\lambda_9$',...
    '$\lambda_{10}$'},...
    'Location','northeast','Interpreter','latex','NumColumns',3)
axis([0 400 1e-15  1e0])
xticks([0:100:400])
yticks([10.^[-15:5:0]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ell_err10_2
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
markers = {'o', '+', '*', 'x', 'square','diamond','^','v','<','>','pentagram','hexagram','-','|','.'};
%
for j=1:10
    semilogy(nv,reerr{3}(:,j),'Marker',markers(j),'LineWidth',1.5)
    hold on; box on
end
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Relative Error','FontSize',14,'Interpreter','latex');
legend({'$\lambda_1$','$\lambda_2$','$\lambda_3$','$\lambda_4$',...
    '$\lambda_5$','$\lambda_6$','$\lambda_7$','$\lambda_8$','$\lambda_9$',...
    '$\lambda_{10}$'},...
    'Location','northeast','Interpreter','latex','NumColumns',3)
axis([0 400 1e-15  1e0])
xticks([0:100:400])
yticks([10.^[-15:5:0]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ell_err10_3
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
markers = {'o', '+', '*', 'x', 'square','diamond','^','v','<','>','pentagram','hexagram','-','|','.'};
%
for j=1:10
    semilogy(nv,reerr{4}(:,j),'Marker',markers(j),'LineWidth',1.5)
    hold on; box on
end
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Relative Error','FontSize',14,'Interpreter','latex');
legend({'$\lambda_1$','$\lambda_2$','$\lambda_3$','$\lambda_4$',...
    '$\lambda_5$','$\lambda_6$','$\lambda_7$','$\lambda_8$','$\lambda_9$',...
    '$\lambda_{10}$'},...
    'Location','northeast','Interpreter','latex','NumColumns',3)
axis([0 400 1e-15  1e0])
xticks([0:100:400])
yticks([10.^[-15:5:0]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ell_err10_4
% 
%%