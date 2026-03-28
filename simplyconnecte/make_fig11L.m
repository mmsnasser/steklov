clc;clear
%%
% plots error; fig 1 plots overall error, fig 2 plots error for each eigenvalue 
resb = [
           20     0.000921       14.219            1
           30    0.0031767       22.414            4
           40    0.0039556       31.432            5
           50    0.0037976       41.356            6
           60     0.003925       51.518            6
           70    0.0057018       61.753            6
           80    0.0046528       72.106            6
           90    0.0050055       82.576            6
          100     0.004906       93.135            6
          110    0.0050561       103.76            6
          120    0.0051806       114.45            6
          130    0.0051506       125.19            6
          140    0.0047221       135.97            6
          150    0.0056615        146.8            6
          160    0.0056119       157.65            6
          170     0.005496       168.54            6
          180    0.0057167       179.46            6
          190    0.0083117        190.4            6
          200    0.0056076       201.37            6
          300     0.010155       311.98            6
          400     0.012358       423.71            6
          500      0.01621       536.13            6
          600     0.024877          649            6
          700     0.032744       762.22            6
          800     0.033108       875.71            6
          900     0.045325        989.4            6
         1000     0.050428       1103.3            6];
%%
resu = [   20    0.0010083       12.311            1
           30    0.0035166       21.927            4
           40    0.0042458       31.797            5
           50    0.0045082       41.516            5
           60    0.0045419       51.489            5
           70    0.0041267       61.724            5
           80    0.0048041       72.105            6
           90    0.0043688       82.579            5
          100    0.0049049       93.135            6
          110    0.0050487       103.76            5
          120    0.0046448       114.45            5
          130    0.0073576       125.19            6
          140    0.0066089       135.97            5
          150    0.0062337        146.8            6
          160    0.0071715       157.65            6
          170    0.0064629       168.54            6
          180    0.0065382       179.46            6
          190    0.0065443        190.4            6
          200    0.0062824       201.37            6
          300    0.0097991       311.98            6
          400     0.012126       423.71            6
          500     0.016158       536.13            6
          600     0.022657          649            6
          700     0.027699       762.22            6
          800     0.029107       875.71            5
          900      0.03723        989.4            6
         1000     0.049297       1103.3            6
         ];
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
loglog(resb(:,1),resb(:,2),'-vb','LineWidth',1.5)
hold on
loglog(resu(:,1),resu(:,2),'-^r','LineWidth',1.5)
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
print -depsc fig_new_ext_time
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
loglog(resb(:,1),resb(:,3),'-vb','LineWidth',1.5)
hold on
loglog(resu(:,1),resu(:,3),'-^r','LineWidth',1.5)
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
print -depsc fig_new_ext_cnd
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
semilogx(resb(:,1),resb(:,4),'-vb','LineWidth',1.5)
hold on
semilogx(resu(:,1),resu(:,4),'-^r','LineWidth',1.5)
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
print -depsc fig_new_ext_itr
% 
%%
