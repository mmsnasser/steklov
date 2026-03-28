clear;clc
% eigenvalues of KD using E
addpath ../bie; addpath ../fmm; 
%%
n =2^10;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
MD = matD(n);
%%
rv = [1;2;2.6;3];
for kj = 1:length(rv)
r = rv(kj);
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
%
et    =  a*(cos(t)+r*1i*sin(t));
etp   =  a*(-sin(t)+r*1i*cos(t));
etpp  =  a*(-cos(t)-r*1i*sin(t));
%
Area  = h*sum(real(et).*imag(etp));
%
A        =  et-alpha;
Ap       =  etp;
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
M        =  -K+tM;
E        =  -inv(eye(n)-N)*M;
%
P  = diag(1./abs(etp));
KD = P*MD*E;
%%
KDI = KD+eye(n);
[V,D] = eigs(KDI,10,'smallestabs');
Dv    = real(diag(D)-1);
%%
[x,y] = meshgrid(linspace(-1,1,250),linspace(-r,r,250));
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
figs(kj) = figure;
tiledlayout(2,4,'TileSpacing','Compact');
% 
for kk=1:8
    kk
    gam  = real(V(:,kk+2)); c = sqrt(h*sum(gam.^2.*abs(etp))); gam = gam./c;
    mun = real(E*gam);
    Fun= gam+i*mun;
    %
    Fzv23 = fcau(et,etp,Fun,zv2);
    uzv23 = real(Fzv23);
    uzv3  = NaN(size(zv));
    uzv3(isfinite(zv)) = uzv23;
    uz3 = reshape(uzv3,size(z));
    % 
%%
nexttile
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
pcolor(real(z),imag(z),uz3)
hold on
vet = et; vet(end+1)=vet(1);
plot(real(vet),imag(vet),'-k',LineWidth=2)
% title("$\lambda$"+kk)
colormap jet
shading interp 
axis equal
axis([min(real(vet))-0.01 max(real(vet))+0.01 min(imag(vet))-0.02 max(imag(vet))+0.02])
axis off
ax=gca; 
set(gca,'LooseInset',get(gca,'TightInset'))
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
% set(gca,'FontSize',14)
drawnow
%% 
end
% print(fig1, 'fig_hypdisk1.pdf', '-dpdf', '-fillpage');
print(figs(kj),'-depsc',sprintf('fig_ex_ell_r%d.eps',kj))
%%
end
%%