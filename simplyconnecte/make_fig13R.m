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
% gam4  = real(V(:,4));
% gam8  = real(V(:,8));
% %% 
% mu4 = real(E*gam4);
% mu8 = real(E*gam8);
% %
% cv4   = gam4+E*mu4; c4 = mean(cv4);
% cv8   = gam8+E*mu8; c8 = mean(cv8);
% % 
% %
% F4= gam4+i*mu4;
% F8= gam8+i*mu4;
%%
[x,y] = meshgrid(linspace(-3,3,250),linspace(-3.5,3.5,250));
xv = real(et);
yv = imag(et);
[in on] = inpolygon(x,y,xv,yv);
x(in)=NaN+i*NaN;
y(in)=NaN+i*NaN;
x(on)=NaN+i*NaN;
y(on)=NaN+i*NaN;
%
z=x+i*y;
zv   = z(:).';
zv2  = zv(isfinite(z));
% 
figure;
tiledlayout(2,4,'TileSpacing','tight','Padding','tight');
% 
for kk=1:8
    gam  = real(V(:,kk+2)); c = sqrt(h*sum(gam.^2.*abs(etp))); gam = gam./c; 
    mun  = real(E*gam);
    Fun  = gam+i*mun;
    cv   = gam+E*mun; cinf = mean(cv);
    %
    Fzv23 = fcau(et,etp,Fun,zv2,n,cinf);
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
% axis equal
axis off
axis ([-3 3 -3.5 3.5])
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
print -depsc fig_kite_8eigf_U