clear;clc
% eigenvalues of KD using E
%%
rv = 0:0.025:0.95; rv = rv(:);
rk = length(rv);
%%
for m=1:rk
    m
    r = rv(m);
    if rv(m)<=.6
        n =2^10;
    else
        n=2^11;
    end
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
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
    I     =  h*sum(sqrt(4*r^2*sin(2*t).^2+(1+r*cos(2*t)).^2));
    a     =  2*pi/I;   
    et    =  a*(1+r*cos(2*t)).*exp(i*t);
    etp   =  a*(-2*r*sin(2*t)+i*(1+r*cos(2*t))).*exp(i*t);
    etpp  =  a*(-5*r*cos(2*t)-1-i*4*r*sin(2*t)).*exp(i*t);
    %
    A        =  et-alpha;
    Ap       =  etp;
    %
%%
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
    plot(rv,xe(k,:),'LineWidth',1.5);
    % plot(rv,xe(k,:),'LineWidth',1.5);
    hold on; box on;
end
legend({'$\lambda_1$','$\lambda_2$','$\lambda_3$','$\lambda_4$',...
    '$\lambda_5$','$\lambda_6$','$\lambda_7$','$\lambda_8$','$\lambda_9$',...
    '$\lambda_{10}$'},...
    'Location','northwest','Interpreter','latex','NumColumns',5)
% set(fig1,'PaperSize',[5  5]);
grid on; axis square
axis([0 0.95 0 7])
xticks([0:0.2:1.0])
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
xlabel('$r$','Interpreter','latex')
ylabel('Eigenvalues')
print -depsc fig_ex_star2_10eig
%%