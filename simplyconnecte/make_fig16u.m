clear;clc
% eigenvalues of KD using E
%%
rv = [1 2 5 10]';
nv = [20:20:400]';
%
xexact = [
         0.999999999999935         0.855836984189535         0.768698829183007         0.747469666437909
         0.999999999999939            1.164010635607          1.25883764854951          1.27303334291789
          1.99999999999996          1.92981202485317           1.8123098096397          1.77436154105265
          1.99999999999997           2.0714075466286          2.17391932481067          2.19304905539063
          2.99999999999996           2.9581957199569          2.82835097834019           2.7775698101104
          2.99999999999998          3.05216718049473          3.21376771973638          3.25193184547303
          3.99999999999996           3.9768183892915          3.84867335440756          3.78819507109412
          3.99999999999999          4.02539858394526          4.16576866101255          4.20766701684323
          4.99999999999995          4.98640358654896          4.86242621841304          4.79363796181362
          4.99999999999997          5.01611374300649          5.17898061757714           5.2401268717324
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
alpha = 0;
%
    I     =  h*sum(sqrt(cos(t).^2+r^2*sin(t).^2));
    a     =  2*pi/I;   
    et    =  a*( cos(t)-r*1i*sin(t));
    etp   =  a*(-sin(t)-r*1i*cos(t));
    etpp  =  a*(-cos(t)+r*1i*sin(t));
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
print -depsc fig_ell_err10_U_1
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
print -depsc fig_ell_err10_U_2
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
print -depsc fig_ell_err10_U_3
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
print -depsc fig_ell_err10_U_4
% 
%%