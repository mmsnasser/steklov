clear;clc
%eigenvalues of KD using E
% Num ex Fig 4, table 1 of The exterior Steklov problem of Euclidean
% domains
addpath ../bie; addpath ../fmm; 
%%
n =2^10;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
et    = 1.5*cos(t)+0.7*cos(2*t)-0.4+i*( 1.5*sin(t)-0.3*cos(t));
etp   =-1.5*sin(t)-1.4*sin(2*t)    +i*( 1.5*cos(t)+0.3*sin(t));
etpp  =-1.5*cos(t)-2.8*cos(2*t)    +i*(-1.5*sin(t)+0.3*cos(t));
Area  = h*sum(real(et).*imag(etp));
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
%%
KDI = KD+eye(n);
[V,D] = eigs(KDI,13,'smallestabs');
%%
[x,y] = meshgrid(linspace(-1.6,1.9,501),linspace(-1.6,1.6,501));
xv = real(et);
yv = imag(et);
[in on] = inpolygon(x,y,xv,yv);
x(~in)=NaN+i*NaN;
y(~in)=NaN+i*NaN;
%
z=x+i*y;
zv   = z(:).';
zv2  = zv(isfinite(z));
% 
figure;
tfig = tiledlayout(2,4);
tfig.TileSpacing = 'none'; % Or 'compact'
tfig.Padding = 'none';     % Or 'compact'
% 
for kk=1:8
    gam  = real(V(:,kk+2)); c = sqrt(h*sum(gam.^2.*abs(etp))); gam = gam./c;
    mun  = real(E*gam);
    Fun  = gam+i*mun;
    %
    Fzv23 = fcau(et,etp,Fun,zv2);
    uzv23 = real(Fzv23);
    uzv3  = NaN(size(zv));
    uzv3(isfinite(zv)) = uzv23;
    uz3 = reshape(uzv3,size(z));
    % 
nexttile

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
pcolor(real(z),imag(z),uz3)
hold on
vet = et; vet(end+1)=vet(1);
plot(real(vet),imag(vet),'-k',LineWidth=2)
colormap jet
shading interp 
axis equal
axis([-1.6 1.9 -1.7 1.7])
axis off
ax=gca; 
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',14)


%% 
end
print -depsc fig_kite_8eigf