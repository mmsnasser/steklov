clear;clc
% eigenvalues of KD using E
%%
% r = [0.2:0.2:5]; r = r(:);
r = linspace(1,10,201); r = r(:);
rk = length(r);
%%
for m=1:rk
    m
    if r(m)<=5
        n =2^10;
    else
        n=2^11;
    end
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
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
    MD = matD(n);
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
figure;
markers = {'o', '+', '*', 'x', 'square','diamond','^','v','<','>','pentagram','hexagram','-','|','.'};
for k=1:10
    % plot(r,xe(k,:),'Marker',markers(k),'LineWidth',1.5);
    plot(r,xe(k,:),'LineWidth',1.5);
    hold on; box on;
end
legend({'$\lambda_1$','$\lambda_2$','$\lambda_3$','$\lambda_4$',...
    '$\lambda_5$','$\lambda_6$','$\lambda_7$','$\lambda_8$','$\lambda_9$',...
    '$\lambda_{10}$'},...
    'Location','northwest','Interpreter','latex','NumColumns',5)
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
axis([1 10 0 6.5])
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$r$','Interpreter','latex')
ylabel('Eigenvalues')
print -depsc fig_ex_ell_10eig_U
%%
figure;
hold on; box on;
plot(r,xe(1,:),'LineWidth',1.5);
plot(r,xa,':r','LineWidth',1.5);
grid on; axis square
axis([1 10 0 3])
legend({'$\lambda_1$','$\frac{1}{a\sqrt{r}}$'},...
    'Location','north','Interpreter','latex')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$r$','Interpreter','latex')
ylabel('Eigenvalues')
print -depsc fig_ex_ell_eig1_U
%%