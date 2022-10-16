%% run this file
% ----------------------------- %
clc; clear;
close all;

trials = 10;

rmse_total = [];
for trial = 1:trials
    tic;
    rmse = particle_filter_GMM(0);
    toc;
    rmse_total = [rmse_total; rmse];
    fprintf('%d iter finished!\n', trial);

end

ARMSE = mean(rmse_total,'all');


% ----------------------------- %
%%
% figure;
% 
% g_iter = [1;5;10];
% 
% plot(g_iter,armse,'-d','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% plot(g_iter,[error error error],'-','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
% hold off;
% legend('distributed','centralised');
% xlabel('consensus iterations');
% ylabel('ARMSE');
%%
% figure;
% 
% g_iter = [1;5;10];
% t = g_iter*24*(1+6+6)/9;
% armse_gmm = [5.081 2.251 2.016 1.73 1.702];
% g_iter1 = [50;100;200];
% t1 = g_iter1*100*2/9;
% plot(t,armse,'-d','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% plot(t1,yd,'-d','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
% hold off;
% legend('Optimal GMM(C=1)','Graph Laplacian(m=100)');
% xlabel('Avg.Scalar per Node per Step');
% ylabel('ARMSE');