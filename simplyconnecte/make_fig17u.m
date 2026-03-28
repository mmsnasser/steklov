clear;clc
%%
nv = [20:10:190,200:100:1000]';
% 
timm =[
     0.018542    0.0025816    0.0021751    0.0018787
     0.049888    0.0072939    0.0063632    0.0080898
     0.009894    0.0086861    0.0096748    0.0074843
    0.0097832    0.0074728    0.0068085    0.0088581
     0.009852    0.0088718    0.0089993    0.0089784
    0.0097161    0.0085002    0.0088911    0.0079794
     0.010523    0.0085778    0.0089346    0.0090788
     0.010688    0.0094331    0.0099871    0.0096397
    0.0088213    0.0087362    0.0094134    0.0092983
    0.0094361    0.0077746     0.007736    0.0097754
    0.0086878     0.009774    0.0075955    0.0094453
      0.01188    0.0085013    0.0084534       0.0101
     0.012987      0.01194     0.014398     0.011807
     0.013561    0.0093931     0.010703     0.010448
      0.01237     0.012872     0.010384      0.01042
     0.011008     0.010637     0.010328     0.010811
     0.014013      0.01087     0.010844    0.0096906
     0.014207     0.013254       0.0125     0.011034
     0.018083     0.011324     0.011684     0.012679
     0.022033     0.022911     0.018354     0.020735
     0.030798     0.024609     0.025606     0.025356
     0.031064     0.047077     0.026716     0.032128
     0.043174     0.098487     0.074049     0.033958
     0.092769     0.043433     0.073647      0.04278
      0.10624     0.090295     0.085621      0.24452
      0.11757      0.13424      0.15469      0.43915
       0.1323      0.30304      0.11611      0.26137
     ];
cnd = [
           10       13.251       24.534       54.467
           15       20.538       37.747       52.896
           20       27.955       53.818       94.107
           25       35.451       69.713       114.61
           30           43        86.04       148.83
           35       50.588       102.53       176.86
           40       58.206       119.18       209.05
           45       65.849       135.94       239.77
           50        73.51       152.81       271.78
           55       81.188       169.75       303.58
           60        88.88       186.76       335.85
           65       96.583       203.83        368.2
           70        104.3       220.95       400.79
           75       112.02       238.12       433.49
           80       119.75       255.33       466.34
           85       127.49       272.57        499.3
           90       135.24       289.84       532.37
           95       142.99       307.14       565.54
          100       150.74       324.47       598.79
          150       228.53       498.82       934.88
          200       306.59       674.43       1275.1
          250        384.8       850.78       1617.8
          300       463.12       1027.6       1962.1
          350       541.51       1204.9       2307.7
          400       619.97       1382.4       2654.2
          450       698.48       1560.1       3001.4
          500       777.03       1738.1       3349.2   
];
% 
itr = [
    1  1  1  1
    2  3  4  5
    4  4  5  5
    6  4  5  5
    6  5  5  5
    6  5  5  5
    6  5  6  5
    6  5  6  6
    6  5  5  5
    7  5  5  5
    6  5  5  5
    6  5  5  5
    6  5  5  5
    7  5  5  6
    7  5  5  5
    7  5  5  5
    7  5  5  5
    7  5  5  5
    7  5  5  5
    7  5  5  5
    7  5  6  6
    7  5  5  6
    7  5  5  5
    7  5  5  5
    7  6  5  5
    7  5  5  5
    7  5  5  6
];
%%
clr = ['m','k','b','r'];
mrk = ['<','>','^','v'];
figure;
box on
for mm=1:4
loglog(nv,timm(:,mm),'color',clr(mm),'marker',mrk(mm),'LineWidth',1.5);
hold on
end
legend({'$r=1$','$r=2$','$r=5$','$r=10$'},...
    'Location','northwest','Interpreter','latex','NumColumns',2)
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('CPU time (sec)','FontSize',14,'Interpreter','latex');
axis([1e1 1e3 1e-3  1e-0])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ell_time_U
%%
%%
figure;
 box on
for mm=1:4
loglog(nv,cnd(:,mm),'color',clr(mm),'marker',mrk(mm),'LineWidth',1.5);
hold on;
end
legend({'$r=1$','$r=2$','$r=5$','$r=10$'},...
    'Location','northwest','Interpreter','latex','NumColumns',2)
% set(fig1,'PaperSize',[5  5]);
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Condition number','FontSize',14,'Interpreter','latex');
axis([1e1 1e3 1e1  1e4])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ell_cnd_U
%%
figure;
 box on
for mm=1:4
semilogx(nv,itr(:,mm),'color',clr(mm),'marker',mrk(mm),'LineWidth',1.5);
hold on;
end
legend({'$r=1$','$r=2$','$r=5$','$r=10$'},...
    'Location','northwest','Interpreter','latex','NumColumns',2)
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Number of iterations','FontSize',14,'Interpreter','latex');
axis([1e1 1e3 1  10])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ell_itr_U
%%