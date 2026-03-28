clear;clc
%eigenvalues of KD using fft
% Table 3 of Numerical studeies of the
% Steklov eigenvalue problem via convformal mappings
addpath ../bie; addpath ../fmm; 
%%
n =2^10;
h=2*pi/n;
t = 0:h:2*pi-h; t=t(:);
%%
F = dftmtx(n);
%%
et   = 8+5*exp(i*t)+0.5*exp(6i*t);
etp  = 5i*exp(i*t)+3i*exp(6i*t);
etpp = -5*exp(1i*t)-18*exp(6i*t);
Area = h*sum(real(et).*imag(etp));
alpha = 8;
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
% print(fig1, 'fig_4air.pdf', '-dpdf', '-fillpage');
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
KDI = KD+eye(n);
[V,D] = eigs(KDI,10,'smallestabs');
%%
[x,y] = meshgrid(linspace(3,14,350),linspace(-5.5,5.5,350));
xv = real(et);
yv = imag(et);
[in on] = inpolygon(x,y,xv,yv);
x(~in)=NaN+i*NaN;
y(~in)=NaN+i*NaN;
%
z=x+i*y;
zv   = z(:).';
zv2  = zv(isfinite(zv));
%% 
figure;
tiledlayout(2,4,'TileSpacing','Compact');
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
% axis equal
axis([2.8 13.8 -6.5 6.5])
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
print -depsc fig_ex3e_8eigf