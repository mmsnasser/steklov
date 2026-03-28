clear;clc
%%
nv = [20:10:190,200:100:1000]';
% 
timm =[
    0.0025315    0.0008497    0.0007663    0.0009515
    0.0074283    0.0029763     0.002675    0.0030214
    0.0036347    0.0027972    0.0028555    0.0033767
     0.004199    0.0031063    0.0030202    0.0036549
    0.0051624     0.003353    0.0031515    0.0044254
    0.0046348    0.0033589    0.0030959    0.0049784
    0.0046104     0.003976    0.0036852     0.010889
    0.0044067    0.0033599    0.0040678    0.0049855
    0.0049085    0.0033814    0.0039819    0.0051748
    0.0046658    0.0035202    0.0048697    0.0058305
    0.0045534    0.0035585    0.0037424    0.0050882
    0.0050441    0.0041228    0.0041979    0.0062779
    0.0059578    0.0045508     0.003859    0.0065706
    0.0079469    0.0049231    0.0043277    0.0065033
     0.006218    0.0073848    0.0074078    0.0059906
    0.0065357    0.0085876    0.0061685    0.0065534
    0.0066922    0.0067154    0.0064633    0.0077185
    0.0068872    0.0059969    0.0061252    0.0076329
    0.0073789     0.006782    0.0057744     0.007889
    0.0098058    0.0079631    0.0098095    0.0094408
     0.013896     0.012751     0.013341     0.013477
     0.020111     0.018253     0.019494     0.024563
     0.024997     0.023503     0.023621      0.02585
     0.033492     0.023978     0.029143     0.030407
     0.035421     0.029352     0.034248     0.036356
     0.044639     0.035912     0.042473     0.049185
     0.061296     0.045942     0.050393     0.053577
     ];
cnd = [
           10       13.246       22.101       28.914
           15       20.538       38.535       69.142
           20       27.955       53.601       84.252
           25       35.451        69.77       120.52
           30           43       86.026       145.62
           35       50.588       102.53       178.63
           40       58.206       119.18       208.13
           45       65.849       135.94       240.25
           50        73.51       152.81       271.54
           55       81.188       169.75        303.7
           60        88.88       186.76       335.79
           65       96.583       203.83       368.23
           70        104.3       220.95       400.77
           75       112.02       238.12        433.5
           80       119.75       255.33       466.33
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
    3  3  4  5
    4  4  5  6
    6  4  5  6
    6  5  5  6
    6  5  5  7
    6  5  6  7
    6  5  6  7
    6  5  6  7
    6  5  6  7
    6  5  5  7
    7  5  6  7
    7  5  5  7
    7  5  5  7
    7  5  6  7
    7  5  5  7
    7  5  6  7
    7  5  6  7
    7  5  6  7
    7  5  6  7
    7  5  6  7
    7  5  5  7
    7  5  5  7
    7  5  6  6
    7  5  6  7
    7  5  6  7
    7  5  6  7
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
axis([1e1 1e3 5e-4  1e-1])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ell_time
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
print -depsc fig_ell_cnd
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
print -depsc fig_ell_itr
%%