clc;clear
%%
% plots error; fig 1 plots overall error, fig 2 plots error for each eigenvalue 
%%
reserr =[
           20     0.063368      0.13875       0.1226      0.16655     0.078755      0.27459      0.43073      0.47751       0.2285      0.67696
           40    0.0030144    0.0085991     0.006448      0.01374     0.034521     0.033888     0.055602     0.074486     0.053969      0.12946
           60   0.00012495   0.00046991   0.00028001   0.00088605    0.0031345    0.0032388    0.0055006    0.0093635     0.024001     0.022357
           80   5.0529e-06   2.2746e-05   1.0629e-05   4.7404e-05   0.00019409   0.00022643   0.00038058    0.0008204    0.0029128    0.0026505
          100   2.0235e-07   1.0354e-06   3.7382e-07   2.3053e-06    1.051e-05   1.3195e-05    2.102e-05   5.5969e-05   0.00023169   0.00022607
          120   8.0627e-09   4.5537e-08   1.2375e-08   1.0636e-07   5.2775e-07   6.9518e-07   1.0218e-06   3.3052e-06   1.5445e-05   1.5664e-05
          140   3.2034e-10   1.9595e-09   3.8291e-10   4.7486e-09   2.5237e-08   3.4435e-08   4.5831e-08    1.788e-07   9.2373e-07   9.5673e-07
          160   1.2724e-11   8.3059e-11   1.0742e-11   2.0737e-10   1.1667e-09   1.6365e-09   1.9406e-09   9.1257e-09   5.1268e-08   5.3808e-08
          180   5.0999e-13   3.4846e-12   2.6493e-13   8.9259e-12   5.2609e-11   7.5522e-11   7.8476e-11   4.4692e-10   2.6958e-09    2.856e-09
          200   3.5189e-14   1.4501e-13   3.1566e-14   3.8348e-13   2.3104e-12   3.4224e-12   3.0429e-12   2.1245e-11   1.3611e-10   1.4524e-10
          220   1.5321e-14   1.3217e-14   1.0334e-14   2.0366e-14   9.2893e-14    1.639e-13   1.1434e-13    9.986e-13   6.6548e-12   7.1465e-12
          240   9.2602e-15   4.7889e-15   2.3111e-14   1.2749e-14   2.6749e-15   1.1425e-14   2.8407e-15   5.6858e-14   3.0849e-13    3.508e-13
          260   1.2964e-14   4.9804e-15   1.9541e-14    7.451e-15   1.7509e-14   1.8675e-14   3.5508e-16   2.0994e-14   3.5442e-15   2.4628e-14
          280   2.0036e-14   7.6622e-15   2.3111e-14   1.2253e-14   1.0943e-14   6.5913e-16   7.9894e-15    1.697e-14   7.2302e-15   8.5357e-15
          300   1.6668e-14   2.8733e-15   1.8789e-14   7.9477e-15   1.7022e-14   2.1092e-14   8.8771e-16   9.6221e-15   1.6729e-14   7.8361e-15
          320   2.1383e-14   6.3213e-15   1.3528e-14   1.4405e-14   2.9181e-15   6.5913e-15   1.7754e-15    1.732e-14   9.4986e-15   1.1614e-14
          340    9.597e-15   1.4175e-14   2.8935e-14   9.9346e-16   3.4044e-15   7.6899e-15   1.2605e-14   8.2225e-15   4.8202e-15   1.5392e-15
          360   1.4985e-14   1.3983e-14    2.048e-14   4.9673e-16   8.0248e-15   5.2731e-15    6.569e-15    9.797e-15   9.3568e-15   9.2354e-15
          380   2.1383e-14   1.5324e-15   3.5887e-14   1.0763e-14   1.3131e-14   4.8336e-15   3.9059e-15    1.557e-14   8.6479e-15   5.7371e-15
          400   1.4985e-14   9.9608e-15   1.7286e-14   1.0266e-14   5.8362e-15   4.1745e-15   5.3262e-16   2.3443e-14   6.6632e-15    9.935e-15
];
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
print -depsc fig_ex4e_err10
% 
%%
%%
