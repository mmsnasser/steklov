clear;clc
%eigenvalues of KD using E
% Num ex Fig 4, table 1 of The exterior Steklov problem of Euclidean
% domains
%%
n =2^10;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
et    = 1.5*cos(t)+0.7*cos(2*t)-0.4+i*( -1.5*sin(t)-0.3*cos(t));
etp   =-1.5*sin(t)-1.4*sin(2*t)    +i*( -1.5*cos(t)+0.3*sin(t));
etpp  =-1.5*cos(t)-2.8*cos(2*t)    +i*(1.5*sin(t)+0.3*cos(t));
Area  = -h*sum(real(et).*imag(etp))
alpha = 0;
%%
fig1=figure;
hold on; box on;
crv = et; crv(end+1)=crv(1);
plot(real(crv),imag(crv),'b','LineWidth',2);
plot(real(alpha),imag(alpha),'sr','MarkerFaceColor','r');
set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',14)
set(ax,'LooseInset',get(ax,'TightInset'))
%%
h        =   2*pi/n;
% 
for k=1:n
    for j=1:n
        if (k==j)
            N(k,j)  =  h*(1/pi)*imag(0.5*(etpp(k)/etp(k)));
            tM(k,j) =  h*(1/pi)*real(0.5*(etpp(k)/etp(k)));
        else
            N(k,j)  =  h*(1/pi)*imag((etp(j))/((et(j)-et(k))));
            tM(k,j) =  h*((1/(2*pi))*cot((t(k)-t(j))/2)...
                        +(1/pi)*real((etp(j))/((et(j)-et(k)))));
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
[V,D] = eigs(KDI,13,'smallestabs');
Dv    = real(diag(D)-1)
%%
% 
figure;
tfig = tiledlayout(2,4);
tfig.TileSpacing = 'none'; % Or 'compact'
tfig.Padding = 'none';     % Or 'compact'
% 
for kk=1:8
    gam  = real(V(:,kk+2)); c = h*sum(gam.^2.*abs(etp)); gam = gam./c;
    % 
nexttile
%
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
etv = et; etv(n+1)=etv(1);
gamv = gam; gamv(n+1)=gam(1);
plot3(real(etv),imag(etv),gamv,'b','LineWidth',1.5)
hold on; box on
plot(real(etv),imag(etv),'-.k','LineWidth',1.5)
for kj=1:8:length(et)
    pnt1 = [real(et(kj)) real(et(kj))];
    pnt2 = [imag(et(kj)) imag(et(kj))];
    pnt3 = [gam(kj)      0];
    if gam(kj)>0
        plot3(pnt1,pnt2,pnt3,'k','LineWidth',0.5);
    else
        plot3(pnt1,pnt2,pnt3,'r','LineWidth',0.5);
    end
end
view([1 1 1])
axis off
ax=gca; 
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',14)
% 
end
%%
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_kite_8eigf_U_b
%%