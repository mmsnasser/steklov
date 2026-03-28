clc;clear
%%
% plots error; fig 1 plots overall error, fig 2 plots error for each eigenvalue 
res3 = [
           20    0.0033545       4.1543            1
           30     0.004662       5.9333            2
           40    0.0050581       7.8471            5
           50    0.0052768       9.9133            6
           60    0.0066083       12.092            7
           70    0.0047551       14.346            7
           80    0.0051434       16.654            7
           90    0.0047155       19.001            7
          100    0.0048329       21.378            7
          110    0.0050687       23.777            7
          120    0.0061574       26.196            8
          130    0.0064815        28.63            8
          140    0.0054146       31.077            7
          150    0.0071635       33.536            8
          160    0.0072357       36.005            8
          170    0.0072236       38.483            8
          180    0.0064615       40.969            7
          190     0.007414       43.462            8
          200    0.0078838       45.962            8
          300     0.012026       71.226            8
          400     0.015095         96.8            8
          500     0.022929       122.56            8
          600     0.024976       148.44            8
          700     0.037387       174.41            8
          800      0.03389       200.44            8
          900     0.045588       226.54            8
         1000     0.056807       252.67            8];
%%
res4 = [
           20    0.0010005       50.003            1
           30    0.0027021       77.732            3
           40    0.0032507       108.58            4
           50     0.003696       138.63            4
           60    0.0038912       169.57            5
           70     0.004153       200.52            4
           80    0.0035556       231.61            5
           90    0.0037452       262.76            5
          100    0.0038876       293.97            5
          110    0.0036863       325.23            5
          120      0.00453       356.53            5
          130    0.0048078       387.86            5
          140    0.0053186       419.22            5
          150    0.0056648       450.61            5
          160    0.0052173       482.02            5
          170    0.0049749       513.45            5
          180    0.0057449       544.91            5
          190    0.0058672       576.38            5
          200    0.0050807       607.87            5
          300    0.0085047       923.41            5
          400     0.011116       1239.8            5
          500      0.01499       1556.6            5
          600      0.02425       1873.8            5
          700     0.026005       2191.2            5
          800     0.030267       2508.8            5
          900     0.033037       2826.5            5
         1000     0.041746       3144.4            5
         ];
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
loglog(res3(:,1),res3(:,2),'-vb','LineWidth',1.5)
hold on
loglog(res4(:,1),res4(:,2),'-^r','LineWidth',1.5)
legend({'$G_1$','$G_2$'},'Location','northwest','Interpreter','latex','NumColumns',3)
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('CPU time (sec)','FontSize',14,'Interpreter','latex');
axis([1e1 1e3 1e-3  1e-1])
xticks([10.^[1:1:3]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ex3_time
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
loglog(res3(:,1),res3(:,3),'-vb','LineWidth',1.5)
hold on
loglog(res4(:,1),res4(:,3),'-^r','LineWidth',1.5)
legend({'$G_1$','$G_2$'},'Location','northwest','Interpreter','latex','NumColumns',3)
xlabel('$n$','FontSize',14,'Interpreter','latex');
axis([1e1 1e3 1e0  1e4])
xticks([10.^[1:1:3]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ex3_cnd
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
semilogx(res3(:,1),res3(:,4),'-vb','LineWidth',1.5)
hold on
semilogx(res4(:,1),res4(:,4),'-^r','LineWidth',1.5)
legend({'$G_1$','$G_2$'},'Location','northwest','Interpreter','latex','NumColumns',3)
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Number of iterations','FontSize',14,'Interpreter','latex');
axis([20 1e3 1  10])
xticks([10.^[1:1:3]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ex3_itr
% 
%%
