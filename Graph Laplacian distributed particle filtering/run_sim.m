% run this file!
%
% -------------------------------------------------------------------------
%
% Description:
% this file reproduces the work of paper 
%  "Rabbat, M., Coates, M., & Blouin, S. (2016, August). Graph Laplacian 
%  distributed particle filtering. In 2016 24th European Signal Processing 
%  Conference (EUSIPCO) (pp. 1493-1497). IEEE."
% -------------------------------------------------------------------------
%
%% centralised
clc;
clear
close all;
% rms = particle_filter_GLap(plt,mode,num_sensors,topo,g_iters,m,GLA,if_subset)
rms = []; % rmse = particle_filter_GLap(1,'multiple',9,'distributed',g_iter,m, 1,0);
iters = 100;
for i=1:iters
    tic;
    rms = [rms; particle_filter_GLap(0,'multiple',9,'centralised',0,0,0,0)];
    toc;
    fprintf('%d finished!\n',i);
end
% rms = particle_filter(1,'multiple',20,'centralised',20,100,0);
% rms = particle_filter_GLap(1,'single');
% rms = rms/iters;
figure;
plot(mean(rms,1),'-');
disp(mean(rms,'all'));
%%
clc;
clear
close all;
% rms = particle_filter_GLap(plt,mode,num_sensors,topo,g_iters,m,GLA,if_subset)
rms = []; % rmse = particle_filter_GLap(1,'multiple',9,'distributed',g_iter,m, 1,0);
iters = 1;
for i=1:iters
    tic;
    rms = [rms; particle_filter_GLap(1,'multiple',9,'distributed',120,500,1,0)];
    toc;
    fprintf('%d finished!\n',i);
end
% rms = particle_filter(1,'multiple',20,'centralised',20,100,0);
% rms = particle_filter_GLap(1,'single');
% rms = rms/iters;
figure;
plot(mean(rms,1),'-');

%%%% results
% GLA-free
% N = 500; rms around 4.0185(50x)
% N = 1000; rms around 1.3928(50x)
% GLA
% N = 2000; m = 50; rms aound 5.3877 (50x)
% N = 2000; m = 500; g_iter = 100; rms aound 1.2730 (20x)
% N = 2000; m = 100; g_iter = 100; rms aound 3.3446 (20x)
% N = 1000; m = 500; rms aound 3.3 (50x)
% %% RMS
% 
% clc; clear;
% close all;
% 
% % hyperparameter setup
% iters = 1;
% T = 30; % time span
% 
% % single-sensor case
% rms_single = zeros(1,T);
% mode = 'single';
% plt = 1;
% for iter = 1:iters
%     rms_single = rms_single + particle_filter(plt,T,mode);
% end
% rms_single = sqrt(rms_single/iters);

% %% multi-sensor case
% % num_sensors = 25;
% rms_multiple_4c0 = zeros(1,T);
% rms_multiple_25c0 = zeros(1,T);
% rms_multiple_25d0 = zeros(1,T);
% rms_multiple_25d1 = zeros(1,T);
% mode = 'multiple';
% m = 100;
% plt = 0;
% for iter = 1:iters
% %    rms_multiple_4c0 = rms_multiple_4c0 + particle_filter(plt,T,mode,4,'centralised',0,m,0);
%     rms_multiple_25c0 = rms_multiple_25c0 + particle_filter(plt,T,mode,25,'centralised',0,m,0);
% %     rms_multiple_25d0 = rms_multiple_25d0 + particle_filter(plt,T,mode,4,'distributed',g_ters,m,0);
% %     rms_multiple_25d1 = rms_multiple_25d1 + particle_filter(plt,T,mode,4,'distributed',g_ters,m,1);
% end
% rms_multiple_4c0 = sqrt(rms_multiple_4c0/iters);
% rms_multiple_25c0 = sqrt(rms_multiple_25c0/iters);
% % rms_multiple_25d0 = sqrt(rms_multiple_25d0/iters);
% % rms_multiple_25d1 = sqrt(rms_multiple_25d1/iters);
% 
% % plot
% figure;
% t = 1:T;
% plot(t,rms_single,'-d','MarkerFaceColor','b','LineWidth',1);
% hold on;
% plot(t,rms_multiple_4c0,'-d','Color','#77AC30','MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1);
% hold on;
% plot(t,rms_multiple_25c0,'-d','Color','#A2142F','MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
% hold on;
% plot(t,rms_multiple_25d0,'-d','Color','#7E2F8E','MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1);
% hold on;
% plot(t,rms_multiple_25d1,'-d','Color','#4DBEEE','MarkerFaceColor','#4DBEEE','MarkerEdgeColor','#4DBEEE','LineWidth',1);
% legend('single-sensor','4-sensor(centralised)','25-sensor(centralised)','25-sensor(distributed)','25-sensor(distributed&approximation)');
% % legend('single-sensor','4-sensor(centralised)','4-sensor(distributed)','4-sensor(distributed&approximation)');
% xlabel('time');
% ylabel('error');
% title('RMS position error');
% hold off;


% %% ARMSE
% clc; clear;
% close all;
% 
% iters = 10;
% T = 30;
% 
% error = zeros(3,3);
% m1 = 100;
% m2 = 500;
% m3 = 1000;
% 
% 
% for g_iter = 500:500:1500
%     rms_m1 = zeros(1,T);
%     rms_m2 = zeros(1,T);
%     rms_m3 = zeros(1,T);
%     for iter = 1:iters
%         % rms = particle_filter_GLap(plt,mode,num_sensors,topo,g_iters,m,GLA)
%         rms_m1 = rms_m1 + particle_filter_GLap(0,'multiple',20,'distributed',g_iter,m1, 1);
%         fprintf('%d g_iter, %d iter, %d m finished!\n', g_iter, iter, m1);
%         rms_m2 = rms_m2 + particle_filter_GLap(0,'multiple',20,'distributed',g_iter,m2, 1);
%         fprintf('%d g_iter, %d iter, %d m finished!\n', g_iter, iter, m2);
%         rms_m3 = rms_m3 + particle_filter_GLap(0,'multiple',20,'distributed',g_iter,m3, 1);
%         fprintf('%d g_iter, %d iter, %d m finished!\n', g_iter, iter, m3);
%         fprintf('%d g_iter, %d iter finished!\n', g_iter, iter);
%     end
%     fprintf('%d g_iter finished!\n', g_iter);
%     rms_m1 = sum(rms_m1)/T/iters;
%     rms_m2 = sum(rms_m2)/T/iters;
%     rms_m3 = sum(rms_m3)/T/iters;
%     error(ceil(g_iter/500),1) = rms_m1;
%     error(ceil(g_iter/500),2) = rms_m2;
%     error(ceil(g_iter/500),3) = rms_m3;
% end
% 
% save ARMSE error;
% 
% % plot
% figure;
% g_iter = 500:500:1500;
% t1 = g_iter*m1*2/20;
% t2 = g_iter*m2*2/20;
% t3 = g_iter*m3*2/20;
% plot(t1,error(:,1),'-d','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% plot(t2,error(:,2),'-d','Color','#A2142F','MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
% hold on;
% plot(t3,error(:,3),'-d','Color','#7E2F8E','MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1);
% hold off;
% legend(['m =',num2str(m1)],['m =',num2str(m2)],['m =',num2str(m3)]);
% xlabel('Avg.Scalar per Node per Step');
% ylabel('ARMSE');
%%
% clc; clear;
% close all;
% 
% iters = 20;
% T = 30;
% m = 100;
% 
% error = zeros(4,1);
% 
% 
% k = 1;
% for g_iter = 10:10:40
%     rms_m = zeros(1,T);
%     for iter = 1:iters
%         % rms = particle_filter_GLap(plt,mode,num_sensors,topo,g_iters,m,GLA)
%         % particle_filter_GLap(1,'multiple',9,'distributed',1000,100, 1,0);
%         rms_m = rms_m + particle_filter_GLap(0,'multiple',4,'distributed',g_iter,m, 1, 0);
%         fprintf('%d g_iter, %d iter finished!\n', g_iter, iter);
%     end
%     rms_m = sum(rms_m)/T/iters;
%     error(k,1) = rms_m;
%     k = k+1;
% end
% 
% save ARMSE error;
% 
% % plot
% figure;
% g_iter = 10:10:40;
% t = (g_iter+8)*m*2/4; % averaging+maximizing
% plot(t,error,'-d','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% xlabel('Avg.Scalar per Node per Step');
% ylabel('ARMSE');
%%
clc; clear;
close all;

T = 30;

g_iter = 300;
m = 500;
rmse_total = [];
for trial = 1:100
rmse = particle_filter_GLap(0,'multiple',9,'distributed',g_iter,m, 1,0);
% rmse = particle_filter_GLap(1,'multiple',9,'centralised',1500,100, 0, 0);
rmse_total = [rmse_total; rmse];
fprintf('%d iter finished!\n', trial);
%rms = particle_filter_GLap(plt,mode,num_sensors,topo,g_iters,m,GLA,if_subset)
end

% gossip
% 9-node 100:100:500
% 100 -> 4; 200 -> -4 ; 500 -> - 30;
%%
% figure;
% m1 = 2000;
% m2 = 1000;
% m3 = 500;
% m4 = 100;
% g_iter0 = [20;50;100;200];
% g_iter = [50;100;200];
% t1 = g_iter0*m1*2/9;
% t2 = g_iter*m2*2/9;
% t3 = g_iter*m3*2/9;
% t4 = g_iter*m4*2/9;
% plot(t1,GLAFalse,'-d','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% plot(t2,m1000,'-d','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
% hold on;
% plot(t3,m500,'-d','Color','#7E2F8E','MarkerSize',4,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1);
% hold on;
% plot(t4,m100,'-d','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1);
% hold off;
% legend(['without approximation.'],['m =',num2str(m2)],['m =',num2str(m3)],['m =',num2str(m4)]);
% xlabel('Avg.Scalar per Node per Step');
% ylabel('ARMSE');
%%
% figure;
% p1 = plot(x1,y1,'LineWidth',1);
% hold on;
% plot(x1(1),y1(1),'-d','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% plot(x1(2),y1(2),'-o','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% plot(x1(3),y1(3),'-s','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% p2 = plot(x2,y2,'Color','#A2142F','LineWidth',1);
% hold on;
% plot(x2(1),y2(1),'-d','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
% hold on;
% plot(x2(2),y2(2),'-o','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
% hold on;
% plot(x2(3),y2(3),'-s','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
% hold on;
% p3 = plot(x3,y3,'Color','#7E2F8E','LineWidth',1);
% hold on;
% plot(x3(1),y3(1),'-d','Color','#7E2F8E','MarkerSize',4,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1);
% hold on;
% plot(x3(2),y3(2),'-o','Color','#7E2F8E','MarkerSize',4,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1);
% hold on;
% plot(x3(3),y3(3),'-s','Color','#7E2F8E','MarkerSize',4,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1);
% hold on;
% p4 = plot(x4,y4,'Color','#77AC30','LineWidth',1);
% hold on;
% plot(x4(1),y4(1),'-x','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1);
% hold on;
% plot(x4(2),y4(2),'-d','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1);
% hold on;
% plot(x4(3),y4(3),'-o','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1);
% hold on;
% plot(x4(4),y4(4),'-s','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1);
% hold off;
% legend([p1,p2,p3,p4],{'m = 100','m = 500','m = 1000','without approximation.'});
% xlabel('Avg.Scalar per Node per Step');
% ylabel('ARMSE');
%%
% see comparison.fig in image/file
figure;
plot(x4,y4,'-d','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
hold on;
plot(x3,y3,'--d','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);
hold on;
plot(x2,y2,'-.d','Color','#7E2F8E','MarkerSize',4,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1);
hold on;
plot(x1,y1,':d','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1);
hold off;
legend('z = 100','z = 500','z = 1000','without approximation.');
xlabel('Avg.Scalar per Node per Step');
ylabel('ARMSE');
%% time complexity(sim)
clc;
clear
close all;
rms = []; % rmse = particle_filter_GLap(1,'multiple',9,'distributed',g_iter,m, 1,0);
T_gldpf = [];
iters = 50;
for i=1:iters
    tic;
    rms = [rms; particle_filter_GLap(0,'multiple',9,'distributed',100,400,1,0)];
    toc;
    T_gldpf(i) = toc;
    fprintf('%d finished!\n',i);
end
% save T_gldpf T_gldpf;