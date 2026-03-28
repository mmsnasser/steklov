clc;clear
%%
%%
reserr =[
           20    0.0013499    0.0013499    0.0061294    0.0061294     0.015362     0.015362    0.0031213    0.0031213     0.031632     0.028994
           40   0.00016573   0.00016573     1.11e-05     1.11e-05    0.0003032    0.0003032   0.00011289   0.00011289   0.00010865   0.00053412
           60   1.3786e-05   1.3786e-05   2.1114e-06   2.1114e-06   2.1528e-05   2.1528e-05   8.8092e-06   8.8092e-06   5.3689e-06   3.0305e-05
           80   1.1999e-06   1.1999e-06   2.4625e-07   2.4625e-07     1.83e-06     1.83e-06    7.902e-07    7.902e-07   3.6267e-07   2.4989e-06
          100   1.1163e-07   1.1163e-07   2.4386e-08   2.4386e-08   1.7027e-07   1.7027e-07   7.4901e-08   7.4901e-08   3.0717e-08   2.2918e-07
          120   1.0989e-08   1.0989e-08   2.4575e-09   2.4575e-09   1.6815e-08   1.6815e-08   7.4674e-09   7.4674e-09    2.882e-09   2.2453e-08
          140   1.1282e-09   1.1282e-09   2.5576e-10   2.5571e-10   1.7312e-09   1.7312e-09   7.7363e-10   7.7361e-10   2.8665e-10   2.2998e-09
          160   1.1955e-10   1.1954e-10   2.7415e-11   2.7353e-11   1.8382e-10   1.8382e-10   8.2531e-11   8.2514e-11   2.9646e-11   2.4328e-10
          180   1.3071e-11   1.3086e-11   3.0609e-12   3.0187e-12   2.0015e-11   2.0016e-11   9.0399e-12   9.0093e-12   3.1309e-12   2.6396e-11
          200   1.5252e-12   1.5151e-12   4.1934e-13   3.6013e-13   2.2571e-12   2.2414e-12   1.0358e-12   1.0126e-12    3.151e-13   2.9457e-12
          220   2.6538e-13   2.6538e-13   1.2129e-13   5.1885e-14   2.7822e-13   2.7766e-13   1.3486e-13   1.3879e-13   2.1165e-14   3.4187e-13
          240   1.3489e-13   9.5989e-14   8.5908e-14   2.1945e-14   7.0201e-14   6.2257e-14   6.8714e-14   4.5759e-14   8.4662e-15   6.1322e-14
          260   1.1732e-13    8.846e-14   8.5908e-14   3.7596e-14   4.8771e-14   3.3438e-14   5.5726e-14   3.0657e-14   2.8044e-14   1.9093e-14
          280   1.2736e-13   7.5913e-14   7.4341e-14   1.3779e-14   4.3599e-14   3.7872e-14   4.8779e-14   2.7334e-14   3.7833e-14   1.4151e-14
          300   4.4544e-14   3.1996e-14   7.0938e-14   1.1058e-14   2.4755e-14   1.8659e-14   5.3914e-14   2.9147e-14   3.3071e-14   2.5607e-14
          320   7.8422e-14   3.3251e-14   6.8897e-14   6.9747e-15    3.547e-14   2.3831e-14   5.7538e-14   2.9147e-14   2.7251e-14   2.3585e-14
          340   7.5913e-14   7.0894e-14   6.2773e-14   1.0377e-14   3.9534e-14   2.3462e-14   4.6967e-14   1.9784e-14   1.1112e-14   3.2795e-14
          360   7.3403e-14   5.0818e-14   5.4607e-14   1.3779e-14   4.1012e-14   3.0113e-14   5.0289e-14   2.8845e-14   7.4079e-15   2.1339e-14
          380   1.0226e-13   9.9753e-14   8.2506e-14   4.5761e-14   2.5125e-14    2.494e-14   5.5726e-14   2.7334e-14   1.0583e-14   3.6838e-14
          400    9.097e-14    6.211e-14   8.2506e-14   1.0377e-14   3.7687e-14   2.6787e-14    3.851e-14   2.4314e-14   1.6403e-14   2.8752e-14 ];
%%
nv = reserr(:,1); 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
markers = {'o', '+', '*', 'x', 'square','diamond','^','v','<','>','pentagram','hexagram','-','|','.'};
for j=1:10
    loglog(nv,reserr(:,j+1),'Marker',markers(j),'LineWidth',1.5)
    hold on; box on
end
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Relative Error','FontSize',14,'Interpreter','latex');
legend({'$\lambda_1$','$\lambda_2$','$\lambda_3$','$\lambda_4$',...
    '$\lambda_5$','$\lambda_6$','$\lambda_7$','$\lambda_8$','$\lambda_9$',...
    '$\lambda_{10}$'},...
    'Location','northeast','Interpreter','latex','NumColumns',3)
axis([20 1e3 1e-16  1e-1])
xticks([10.^[1:1:3]])
yticks([10.^[-18:4:-1]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(gca,'LooseInset',get(gca,'TightInset'))
%
% 
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
markers = {'o', '+', '*', 'x', 'square','diamond','^','v','<','>','pentagram','hexagram','-','|','.'};
for j=1:10
    semilogy(nv,reserr(:,j+1),'Marker',markers(j),'LineWidth',1.5)
    hold on; box on
end
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Relative Error','FontSize',14,'Interpreter','latex');
legend({'$\lambda_1$','$\lambda_2$','$\lambda_3$','$\lambda_4$',...
    '$\lambda_5$','$\lambda_6$','$\lambda_7$','$\lambda_8$','$\lambda_9$',...
    '$\lambda_{10}$'},...
    'Location','northeast','Interpreter','latex','NumColumns',3)
axis([0 400 1e-15  1e0])
xticks([0:100:400])
yticks([10.^[-15:5:0]])
grid on; 
axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -depsc fig_ex3e_err10
% 
%%