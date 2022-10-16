%% centralised and single-sensor(test GPR)
clear; clc; close all;
network = 'distributed';
N = 50;
plt = 1;
rmse = particle_filter_GP(N, network, 1, plt);
%% Solution 1.b(run sim)
clear;clc;close all;
trials = 100;
N = [150 500 1000]; % particle size
network = 'centralised';
iters = 0;
armse = [];
for n = 1:length(N)
    rms = [];
    for i=1:trials
        rms = [rms; particle_filter_GP(N(n), network, iters, 0)];
        fprintf('%d finished!\n',i);
    end
    armse = [armse; mean(rms,1)];
end
%% 1.b plot
figure;
plot(armse(1,:)); hold on;
plot(armse(2,:)); hold on;
plot(armse(3,:)); hold off;
legend('M = 150','M = 500','M = 1000');
%% Solution 1.d(run sim)
clear;clc;close all;
trials = 100;
N = [50 100 200]; % particle size
N = 50;
network = 'distributed';
iters = 2;
armse = zeros(length(N),length(iters));
for n = 1:length(N)
    for iter = 1:length(iters)
        rms = [];
        for i=1:trials
            tic;
            rms = [rms; particle_filter_GP(N(n), network, iters(iter), 0)];
            toc;
            fprintf('%d finished!\n',i);
        end
        fprintf('%d n %d iter finished!\n',N(n), iters(iter));
        armse(n, iter) = mean(rms,'all');
    end
end
%%
clear;clc;close all;
trials = 500;
N = 50; % particle size
network = 'distributed';
iters = 0;
armse = [];
for iter = 1:length(iters)
    rms = [];
    for i=1:trials
        rms = [rms; particle_filter_GP(N, network, iters(iter), 0)];
        fprintf('%d finished!\n',i);
    end
    clc; fprintf('%d iter finished!\n',iters(iter));
    armse = [armse; mean(rms,'all')];
end
%% solution 1(plot)
% 1.a & 1.b
x = [50 100 150 200 500]; % particle
y1 = [38.6158 18.2907 12.2316 8.6153 3.8017]; % standard resampling
y2 = [15.5324 6.1984 3.3119 2.7141 2.3]; % GP-based
figure;
plot(x,y1,'-s','Color','b','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);hold on;
plot(x,y2,'--s','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);hold off;
xlabel('particle size'); ylabel('ARMSE');
legend('Standard Resampling', 'GP-enhanced Resampling');

%% 1.d
x = [1 2 3 4 5]; % iterations
x = [0 1 2 3 4 5]; % iterations
y50 = [21.3530 16.9240 11.63 19.6462 26.9929];
y50 = [18.7949 17.0691 16.0998 20.1375 22.6353 30.8394];
y100 = [9.7536 6.5854 6.7392 7.6446 7.1917];
y100 = [8.5111 7.1518 6.2969 7.3154 8.0626 9.2737];
y200 = [5.9661 4.2010 3.3832 3.4796 3.4806];
y200 = [5.8113 4.3797 3.3431 3.2869 3.4010 3.8220];
figure;
plot(x,y50,':s','Color','b','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);hold on;
plot(x,y100,'--s','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1);hold on;
plot(x,y200,'-s','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1);hold off;
xlabel('Iterations'); ylabel('ARMSE');
legend('M = 50', 'M = 100', 'M = 200');
% figure;
% plot(x,y50,'Linewidth',1);
% xlabel('Iterations'); ylabel('ARMSE');
% legend('M = 50');
% figure;
% plot(x,y100,'Linewidth',1);hold on;
% plot(x,y200,'Linewidth',1);hold off;
% xlabel('Iterations'); ylabel('ARMSE');
% legend('M = 100', 'M = 200');
%% Solution 2(sim)
% 2.b
close all; clc; clear;
trials = 50;
N_sensor = 9;
N = 50; % particle size
phi = 1;
M = round(N*phi*N_sensor); % particle size in unit center
iters_p = [80 100 150 200]; iters_w = 20;
% network = 'centralised';
network = 'distributed';
T = []; armse = zeros(1,length(iters_p));
for iter_p = 1:length(iters_p)
    rms = [];
    for i=1:trials
        tic;
        rms = [rms; particle_filter_GP2(N, M, iters_p(iter_p), iters_w, [], network, 'b', 0)];
        toc;
        T(i) = toc;
        fprintf('%d finished!\n',i);
    end
    armse(1,iter_p) = mean(rms,'all');
end

%% fixed iters_w
close all; clc; clear;
trials = 20;
N_sensor = 9;
N = 200; % particle size
phi = [0.07 0.1, 0.2, 0.3, 0.5];
M = round(N*phi*N_sensor); % particle size in unit center
iters_p = [80, 100, 120, 150]; iters_w = [10, 20, 50, 100];
% network = 'centralised';
network = 'distributed';
armse = zeros(length(M), length(iters_p)); 
for M_id = 1:1%length(M)
    for iters_p_id = 2:3
        rms = [];
        for i=1:trials
            rms = [rms; particle_filter_GP2(N, M(M_id), iters_p(iters_p_id), iters_w(2), network, 'b', 0)];
            fprintf('%d finished!\n',i);
        end
        armse(M_id, iters_p_id) = mean(rms,'all');
    end
end

% save armse_fixed_itersw armse;
%% plot
load armse_fixed_itersw;
iters_p = [80, 100, 150, 200];
figure;
x = iters_p;
% for i = 1:size(armse,1)
%     plot(x, armse(i,:),'-s', 'Linewidth',1,'MarkerFaceColor',colors(ID(i),:),'MarkerEdgeColor',colors(ID(i),:),'MarkerSize',4); hold on;
% end
figure; 
plot(x, armse(1,:),'-s','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
plot(x, armse(2,:),'-s','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
plot(x, armse(3,:),'-s','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold on;
plot(x, armse(4,:),'-s','Linewidth',1,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerSize',4); hold off;

hold off;
legend('\phi = 0.1','\phi = 0.2','\phi = 0.3','\phi = 0.5');
xlabel('gossip iterations(particles)'); ylabel('ARMSE');

%% fixed iters_p
close all; clc; clear;
trials = 100;
N_sensor = 9;
N = 100; % particle size
phi = [0.1, 0.2, 0.3, 0.5];
M = N*phi*N_sensor; % particle size in unit center
iters_p = [80, 100, 150, 200]; iters_w = [0, 1, 5, 10, 20, 50];
% network = 'centralised';
network = 'distributed';
armse = zeros(length(M), length(iters_w)); 
for M_id = 2:2%length(M)
    for iters_w_id = 1:length(iters_w)
        rms = [];
        for i=1:trials
            rms = [rms; particle_filter_GP2(N, M(M_id), iters_p(end), iters_w(iters_w_id), network, 'b', 0)];
            fprintf('%d finished!\n',i);
        end
        armse(M_id, iters_w_id) = mean(rms,'all');
    end
end
%% save armse_fixed_itersp armse;
load armse_fixed_itersp;
iters_w = [0, 1, 5, 10, 20, 50];
figure;
x = iters_w;
% for i = 1:size(armse1,1)
%     plot(x, armse1(i,:), 'Linewidth',1); hold on;
% end
% for i = 1:size(armse1,1)
%     semilogy(x, armse1(i,:), 'Linewidth',1); hold on;
% end

figure; semilogy(x, armse1(1,:),'-s','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
semilogy(x, armse1(2,:),'-s','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
semilogy(x, armse1(3,:),'-s','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold on;
semilogy(x, armse1(4,:),'-s','Linewidth',1,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerSize',4); hold off;

hold off;
grid on;
legend('\phi = 0.1','\phi = 0.2','\phi = 0.3','\phi = 0.5');
xlabel('gossip iterations(weights)'); ylabel('ARMSE(log scale)');
%%
close all; clc; clear;
trials = 50;
N_sensor = 9;
N = 50; % particle size3
phi = 1;
M = N*phi*N_sensor; % particle size in unit center
iters_p = 150; iters_w = 50;
% network = 'centralised';
network = 'distributed';
rms = [];
for i=1:trials
    tic;
    rms = [rms; particle_filter_GP2(N, M, iters_p, iters_w, [], network, 'b', 0)];
    toc;
    % particle_filter_GP2(N, M, iters_p, iters_w, iters, network, distributed, plt)
    fprintf('%d finished!\n',i);
end
%% 2.c
clear;clc;close all;
trials = 5;
N_sensor = 9;
iters = [0 2 3 4 5];
N = 50; % particle size
phi = 0.5;
M = round(N*phi*N_sensor); 
rms = [];
for i=1:trials
    tic
    rms = [rms; particle_filter_GP2(N, M, [], [], iters, 'distributed', 'c', 1)];
    toc;
    fprintf('%d finished!\n',i);
end
%% 2.c
clear;clc;close all;
trials = 50;
N_sensor = 9;
iters = [0 1 2 3 4 5];
N = 50; % particle size
phi = 1;
M = round(N*phi*N_sensor); % particle size in unit center
armse = zeros(length(N),length(iters));
for n = 1:length(N)
    for iter = 1:length(iters)
        rms = [];
        for i=1:trials
            tic;
            rms = [rms; particle_filter_GP2(N(n), M(n), [], [], iters(iter), 'distributed', 'c', 0)];
            toc;
            fprintf('%d finished!\n',i);
        end
        armse(n,iter) = mean(rms,'all');
    end
end


%% Solution 2(plot)
% 2.a
phi = [0.1,0.2,0.3,0.5];
N50 = [20.6, 3.6037, 1.9189, 1.6586];
N70 = [7.927, 1.7573, 1.5682, 1.4830];
N100 = [2.0239, 1.5593, 1.4288, 1.4009];
N200 = [1.4954, 1.3186, 1.3091, 1.2883];
% figure; plot([0.1,0.2,0.3,0.5], N50,'-','Linewidth',1,'MarkerFaceColor','auto'); hold on;
% plot(phi, N70,'-','Linewidth',1); hold on;
% plot(phi, N100,'-','Linewidth',1); hold on;
% plot(phi, N200,'-','Linewidth',1); hold off;
figure; semilogy(phi, N50,'-s','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
semilogy(phi, N70,'-s','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
semilogy(phi, N100,'-s','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold on;
semilogy(phi, N200,'-s','Linewidth',1,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerSize',4); hold off;
% ylim([0 100]);
grid on;
legend('M = 50','M = 70','M = 100', 'M = 200'); xlabel('\phi'); ylabel('ARMSE(log scale)');

%% 2.c
x = [0 1 2 3 5]; % iterations
y100 = [3.31  2.11  1.9537  1.8992  1.8474];
y200 = [2.7716 1.96155 1.78785 1.77775 1.74375];
y300 = [2.6845 1.9198 1.7529 1.6957 1.6966];
figure;
plot(x, y100,'-s','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
plot(x, y200,'-s','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
plot(x, y300,'-s','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold off;
xlabel('Iterations'); ylabel('ARMSE');
legend('M = 100', 'M = 200', 'M = 300');
xlabel('Iterations'); ylabel('ARMSE');


%% GPGA(sim)
% apply the GA in distributed sense
clear;clc;close all;
trials = 5;
protocol = 'broadcast';
crossover = 'MFX';
protocol = 'gossip';
N = 20; % particle size
iters = 100;
rms = [];
for i=1:trials
%     tic;
    rms = [rms; particle_filter_GPGA(N, iters, protocol, crossover, 1)];
%     toc;
    fprintf('%d finished!\n',i);
end

%%
clear;clc;close all;
trials = 50;
protocol = 'broadcast';
crossover = 'MFX';
N = 5; % particle size
iters = [5 10 20 30 40 50];
iters = 80;
armse_20 = [];
for iter = 1:length(iters)
    rms = [];
    for i=1:trials
        tic;
        rms = [rms; particle_filter_GPGA(N, iters(iter), protocol, crossover, 0)];
        toc;
        fprintf('%d finished!\n',i);
    end
    armse_20 = [armse_20 mean(rms,'all')];
end

%% GPGA(plot)
% gossip based
iter = [100 150 200 300 500];
armse_10 = [12.2051 5.9199 2.2774 1.8911 1.4348];
armse_20 = [5.7190 2.2561 1.9679 1.6327 1.2815];
armse_50 = [2.9412 2.3751 1.8675 1.5868 1.2678];
figure;
plot(iter, armse_10,'-^','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
plot(iter, armse_20,'-^','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
plot(iter, armse_50,'-^','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold off;
% plot(iter, armse_10, 'Linewidth', 1); hold on;
% plot(iter, armse_20, 'Linewidth', 1); hold on;
% plot(iter, armse_50, 'Linewidth', 1); hold off;
xlabel('Gossip Iterations'); ylabel('ARMSE');
legend('M = 10','M = 20','M = 50');
% broadcast based
iter = [5 10 20 30 40 50];
armse_5 = [9.7765 3.6608 3.1044 2.2360 1.9467 1.7344];
armse_10 = [2.7846 2.2366 2.0315 1.8323 1.7326 1.6941];
armse_20 = [2.3774 2.0791 1.8945 1.7837 1.6868 1.6184];
figure;
plot(iter, armse_5,'-^','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
plot(iter, armse_10,'-^','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
plot(iter, armse_20,'-^','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold off;
% plot(iter, armse_5, 'Linewidth', 1); hold on;
% plot(iter, armse_10, 'Linewidth', 1); hold on;
% plot(iter, armse_20, 'Linewidth', 1); hold off;
xlabel('Broadcast Iterations'); ylabel('ARMSE');
legend('M = 5','M = 10','M = 20');
%% GPFA(sim)
% apply the FA in distributed sense
clear;clc;close all;
trials = 100;
protocol = 'gossip';
% protocol = 'broadcast';
% N = [5 10 20 30 50]; % particle size
N = 10;
iters_FA = 1;
%iters_com = [50 80 100 150 200];
iters_com = 50;
armse = [];
for iter = 1:length(iters_com)
    for i=1:trials
        tic;
        rms = particle_filter_GPFA(N, iters_FA, iters_com(iter), protocol, 0);
        toc;
        fprintf('%d finished!\n',i);
        armse = [armse mean(rms,'all')];
    end
end

%%
% apply the FA in distributed sense
clear;clc;close all;
trials = 500;
protocol = 'gossip';
% protocol = 'broadcast';
% N = [5 10 20 30]; % particle size
N = 5;
iters_FA = 1;
iters_com = [40 50 70 100 150];
% iters_com = [40 150];
% iters_com = 40;
armse = [];
for iter_com = 1:length(iters_com)
    for n = 1:length(N)
        rms = [];
        for i=1:trials
    %         tic;
            rms = [rms; particle_filter_GPFA(N(n), iters_FA, iters_com(iter_com), protocol, 0)];
    %         toc;
            fprintf('%d finished!\n',i);
        end
        armse(iter_com, n) = mean(rms,'all');
    end
end

%% GPFA(plot)
% iters = [50 80 100 150 200];
% y5 = [4.6607 2.0261 1.8090 1.8664 1.9282];
% y10 = [1.8217 1.6982 1.7342 1.7869 1.8490];
% y20 = [1.7351 1.6763 1.6835 1.7527 1.8167];
% y30 = [1.7016 1.6508 1.6829 1.7470 1.7984];
% y50 = [1.6781 1.6374 1.6595 1.7288 1.7874];
% figure;
% % plot(iters, y5, 'Linewidth', 1); hold on;
% plot(iters, y10, 'Linewidth', 1); hold on;
% plot(iters, y20, 'Linewidth', 1); hold on;
% plot(iters, y30, 'Linewidth', 1); hold on;
% plot(iters, y50, 'Linewidth', 1); hold off;
% % legend('M = 5','M = 10','M = 20','M = 30','M = 50');
% legend('M = 10','M = 20','M = 30','M = 50');
% xlabel('Gossip iterations'); ylabel('ARMSE');
close all;clear;clc;
iters = [40 50 70 100 150];
y5 = [11.0804 2.9116 1.7868 1.6359 1.6291];
y10 = [2.3228 1.7453 1.5625 1.4977 1.4982];
y20 = [1.6715 1.5401 1.4640 1.4330 1.4500];
y30 = [1.6763 1.5178 1.4258 1.4262 1.4433];
figure;
% semilogy(iters, y5, 'Linewidth', 1); hold on;
% semilogy(iters, y10, 'Linewidth', 1); hold on;
% semilogy(iters, y20, 'Linewidth', 1); hold on;
% semilogy(iters, y30, 'Linewidth', 1); hold off;
semilogy(iters, y5,'-h','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
semilogy(iters, y10,'-h','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
semilogy(iters, y20,'-h','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold on;
semilogy(iters, y30,'-h','Linewidth',1,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerSize',4); hold off;
xlim([40 150]);
grid on;
xlabel('Gossip iterations'); ylabel('ARMSE');
legend('M = 5','M = 10','M = 20','M = 30');
y5 = [11.0804 2.9116 1.7868 1.6359 1.6291];
y10 = [2.3228 1.7453 1.5625 1.4977 1.4982];
y20 = [1.6715 1.5401 1.4640 1.4330 1.4500];
y30 = [1.6763 1.5178 1.4258 1.4262 1.4433];
figure;
% plot(iters, y5, 'Linewidth', 1); hold on;
plot(iters, y10,'-h','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
plot(iters, y20,'-h','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
plot(iters, y30,'-h','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold off;
% plot(iters, y10, 'Linewidth', 1); hold on;
% plot(iters, y20, 'Linewidth', 1); hold on;
% plot(iters, y30, 'Linewidth', 1); hold off;
xlim([40 150]);
xlabel('Gossip iterations'); ylabel('ARMSE');
legend('M = 10','M = 20','M = 30');

%% Comparison
%clear;clc;close all;
iter1 = [50 100 200 300];
com11 = iter1*100*2/9;
com12 = iter1*500*2/9;
com1 = iter1*200*2/9;
y1 = [ 2.0643 1.8959];
y11 = [7.852 3.737 3.18];
y12 = [6.176 2.272 1.343 1.7333];
iter2 = [2 4 6 8 10];
com2 = 28*24*iter2/9;
y2 = [5.081 2.251 2.016 1.73 1.702];
iter3 = [1 2 3 4 5];
com3 = (24+iter3*80*97)/9;
y3 = [5.9661 4.2010 3.3832 3.4796 3.4806];
iter4 = [80 100 150 200]; % Lp
com4 = (7200+iter4*240*9)/9;%M=200,Lw=20, phi=0.1
y4 = [3.8320 2.3811 1.6900 1.6708];
com4 = (18000+iter4*600*9)/9;%M=50,Lw=20, phi=1
y4 = [3.8314 2.3094 1.6706 1.6603];
% iter4 = [80 100 120 150];
% com4 = (200*0.07*9*2*20+iter4*12*14*9)/9; %M=200,Lw=20, phi=0.07
% y4 = [5.5413 2.5198 1.8804 1.7717];
iter5 = [1 2 3 4 6];
iter5 = [1 2 3 4 5 6];
N = 50;
com5 = iter5*(N*0.5*6*24)/9;
com5 = iter5*(N*1*6*24)/9; % N = 50; phi = 1
y5 = [2.7716 1.96155 1.78785 1.77775 1.74375]; % N = 200;
y5 = [3.1709 2.2954 2.0683 1.9903 1.9785 1.9699]; % N = 50;
y5 = [3.8654 2.2340 1.9622 1.8952 1.8473 1.8612]; % N = 50; phi = 1
iter6 = [100 150 200 300 500];
N = 20;
% iter6 = [300 350 400 500 600 1000];
com6 = 2*iter6*N*6/9;
% y6 = [12.2051 5.9199 2.2774 1.8911 1.4348]; % GPGA(M = 10);
% y6 = [2.8121 2.1037 1.9202 1.7158 1.5740 1.4144]; % GPGA(M = 5);
y6 = [5.7190 2.2561 1.9679 1.6327 1.2815]; % (M = 20);
iter7 = [8 16 24 40];
iter7 = [5 10 20 30 40 50];
M = 10;
com7 = 2*12*iter7*M*6/9;
y7 = [4.4063 3.9332 3.6181 2.9790]; % GPGA(broadcast)(CMX)(M = 10);
y7 = [9.7765 3.6608 3.1044 2.2360 1.9467 1.7344];% GPGA(broadcast)(MFX)(M = 5);
y7 = [2.7846 2.2366 2.0315 1.8323 1.7326 1.6941];% GPGA(broadcast)(MFX)(M = 10);
%iter8 = 50;
iter8 = 80;
N = [3 5 10 20 30 50];
com8 = iter8*N*(2*6+4)/9;
%y8 = [4.6607 1.8217 1.7351 1.7016 1.6781];
y8 = [2.4229 2.0261 1.6982 1.6763 1.6508 1.6374];
iters8 = [40 50 70 100 150];
N = 10;
com8 = iters8*N*(2*6+4)/9;
y8 = [2.3228 1.7453 1.5625 1.4977 1.4982];
% com8(1) = 70*5*(2*6+4)/9;
% com8 = [50*5*(2*6+4)/9 com8];
% y8 = [2.9116 1.7868 1.6359 1.5625 1.4977 1.4982];
%%
figure;
semilogx(com2,y2,'-o','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1); hold on;
semilogx(com11,y11,'--d','Color','#0072BD','MarkerSize',4,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1); hold on;
semilogx(com12,y12,':d','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1); hold off;
xlim([10^2 10^5]);
grid on;
xlabel('Avg.scalar per Node per Step'); ylabel('ARMSE');
% legend('GLDPF(m=500)','GMMDPF(C=1)','GPDPF(solution1.d)(M=200)','GPDPF(solution2.b)(M=200,\phi=0.07)','GPDPF(solution2.c)(M=200,\phi=0.5)');
legend('GMMDPF(C=1)','GLDPF(z=100)','GLDPF(z=500)');
%% conference plot
figure;
semilogx([10^2,10^5],[1.2,1.2],'--'); hold on;
semilogx(com2,y2,'-.s','Color','#0072BD','MarkerSize',4,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1); hold on;
% semilogx(com3,y3,'--s','Color','#7E2F8E','MarkerSize',4,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1); hold on;
semilogx(com8,y8,'--+','Color','#D95319','MarkerSize',4,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','LineWidth',1); hold on;
semilogx(com5,y5,'->','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1); hold on;
% semilogx(com7,y7,'-.^','Color','#EDB120','MarkerSize',4,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','LineWidth',1); hold on;

% semilogx(com6,y6,'--^','Color','#4DBEEE','MarkerSize',4,'MarkerFaceColor','#4DBEEE','MarkerEdgeColor','#4DBEEE','LineWidth',1); hold on;
semilogx(com12,y12,':x','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1); hold on;
% semilogx(com4,y4,'-.s','Color','#FF00FF','MarkerSize',4,'MarkerFaceColor','#FF00FF','MarkerEdgeColor','#FF00FF','LineWidth',1); hold on;
xlim([10^2 1*10^5]);
grid on;
xlabel('Avg.scalar per Node per Step'); ylabel('ARMSE');
% legend('GLDPF(m=500)','GMMDPF(C=1)','GPDPF(solution1.d)(M=200)','GPDPF(solution2.b)(M=200,\phi=0.07)','GPDPF(solution2.c)(M=200,\phi=0.5)');
% legend('GM-DPF(C=1)','GPDPF(Algorithm.8)(M=200)','GPDPF(diffusion)(M=50,\phi=1)','GADPF(M = 10)(broadcast)','FADPF(M = 10)(gossip)','GADPF(M = 20)(gossip)','GLDPF(z=500)','GPDPF(consensus)(M=50,\phi=1)');
legend('Centralised','GM-DPF','FA-DPF','GP-DPF','GL-DPF');
legend('Location','northwest');
%%
figure;
semilogx(com2,y2,':o','Color','#0072BD','MarkerSize',4,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineWidth',1); hold on;
semilogx(com3,y3,'--s','Color','#7E2F8E','MarkerSize',4,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineWidth',1); hold on;
semilogx(com5,y5,'-s','Color','#77AC30','MarkerSize',4,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineWidth',1); hold on;
semilogx(com7,y7,'-.^','Color','#EDB120','MarkerSize',4,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','LineWidth',1); hold on;
semilogx(com8,y8,'-h','Color','#D95319','MarkerSize',4,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','LineWidth',1); hold on;
semilogx(com6,y6,'--^','Color','#4DBEEE','MarkerSize',4,'MarkerFaceColor','#4DBEEE','MarkerEdgeColor','#4DBEEE','LineWidth',1); hold on;
semilogx(com12,y12,':d','Color','#A2142F','MarkerSize',4,'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','LineWidth',1); hold on;
semilogx(com4,y4,'-.s','Color','#FF00FF','MarkerSize',4,'MarkerFaceColor','#FF00FF','MarkerEdgeColor','#FF00FF','LineWidth',1); hold on;
xlim([10^2 3*10^5]);
grid on;
xlabel('Avg.scalar per Node per Step'); ylabel('ARMSE');
% legend('GLDPF(m=500)','GMMDPF(C=1)','GPDPF(solution1.d)(M=200)','GPDPF(solution2.b)(M=200,\phi=0.07)','GPDPF(solution2.c)(M=200,\phi=0.5)');
legend('GMMDPF(C=1)','GPDPF(Algorithm.8)(M=200)','GPDPF(diffusion)(M=50,\phi=1)','GADPF(M = 10)(broadcast)','FADPF(M = 10)(gossip)','GADPF(M = 20)(gossip)','GLDPF(z=500)','GPDPF(consensus)(M=50,\phi=1)');
%% time complexity(sim)(50x)
% GPDPF(consensus)
close all; clc; clear;
trials = 50;
N_sensor = 9;
N = 50; % particle size
phi = 0.6;
M = round(N*phi*N_sensor); % particle size in unit center
iters_p = 150; iters_w = 10;
network = 'distributed';
rms = [];
T_gpdpf_c = [];
for i=1:trials
    tic;
    rms = [rms; particle_filter_GP2(N, M, iters_p, iters_w, [], network, 'b', 0)];
    toc;
    T_gpdpf_c(i) = toc;
    fprintf('%d finished!\n',i);
end
armse = mean(rms,'all');
%save T_gpdpf_c T_gpdpf_c;
%% GPDPF(diffusion)
close all; clc; clear;
trials = 50;
N_sensor = 9;
N = 40; % particle size
phi = 0.5;
M = round(N*phi*N_sensor); % particle size in unit center
iters = 5;
rms = [];
T_gpdpf_d = [];
for i=1:trials
    tic;
    rms = [rms; particle_filter_GP2(N, M, [], [], iters, 'distributed', 'c', 0)];
    toc;
    T_gpdpf_d(i) = toc;
    fprintf('%d finished!\n',i);
end
armse = mean(rms,'all');
%save T_gpdpf_d T_gpdpf_d;
%% GADPF(gossip)
close all; clc; clear;
trials = 50;
N_sensor = 9;
N = 10; % particle size
iters = 250;
rms = [];
T_gadpf_g = [];
for i=1:trials
    tic;
    rms = [rms; particle_filter_GPGA(N, iters, 'gossip', [], 0)];
    toc;
    T_gadpf_g(i) = toc;
    fprintf('%d finished!\n',i);
end
armse = mean(rms,'all');
save T_gadpf_g T_gadpf_g;
%% GADPF(broadcast)
close all; clc; clear;
trials = 50;
N_sensor = 9;
N = 10; % particle size
iters = 20;
rms = [];
T_gadpf_b = [];
for i=1:trials
    tic;
    rms = [rms; particle_filter_GPGA(N, iters, 'broadcast', 'MFX', 0)];
    toc;
    T_gadpf_b(i) = toc;
    fprintf('%d finished!\n',i);
end
armse = mean(rms,'all');
save T_gadpf_b T_gadpf_b;
%% FADPF(gossip)
close all; clc; clear;
trials = 50;
N_sensor = 9;
N = 10; % particle size
rms = [];
T_fadpf_g = [];
for i=1:trials
    tic;
    rms = [rms; particle_filter_GPFA(N, 1, 42, 'gossip', 0)];
    toc;
    T_fadpf_g(i) = toc;
    fprintf('%d finished!\n',i);
end
armse = mean(rms,'all');
save T_fadpf_g T_fadpf_g;
%% time complexity(plot)
clear;
load T_gldpf;
load T_gmmdpf;
load T_gpdpf_c;
load T_gpdpf_d;
load T_gadpf_g;
load T_gadpf_b;
load T_fadpf_g;
figure;
boxplot([T_gldpf', T_gmmdpf', T_gpdpf_c', T_gpdpf_d', T_gadpf_g',T_gadpf_b', T_fadpf_g'],'Labels',{'GL-DPF','GM-DPF','GP-DPF(consensus)','GP-DPF(diffusion)','GA-DPF(gossip)','GA-DPF(broadcast)','FA-DPF(gossip)'});
ylabel('Time/s'); xlabel('Algorithm');
figure;
boxplot([T_gpdpf_c', T_gpdpf_d', T_gadpf_g',T_gadpf_b', T_fadpf_g'],'Labels',{'GP-DPF(consensus)','GP-DPF(diffusion)','GA-DPF(gossip)','GA-DPF(broadcast)','FA-DPF(gossip)'});
ylabel('Time/s'); xlabel('Algorithm');
%% conference plot
figure;
boxplot([T_gldpf', T_gmmdpf', T_gpdpf_d', T_fadpf_g'],'Labels',{'GL-DPF','GM-DPF','GP-DPF','FA-DPF'});
ylabel('Time/s'); xlabel('Algorithm');
%% ARMSE as a function of particle M(sim)
clear;clc;close all;
trials = 100;
N_sensor = 9;
iters_p = 500; iters_w = 50;
iters_c = 5;
iters_ga_g = 500;
iters_ga_b = 50;
iters_fa_g = 100;
N = [10 20 30 50]; % particle size
N = 10;
phi = 1;
M = round(N*phi*N_sensor); % particle size in unit center
armse = zeros(length(N),5);
for n = 1:length(N)
    rms1 = []; rms2 = []; rms3 = []; rms4 = []; rms5= [];
    for i=1:trials
%         rms1 = [rms1; particle_filter_GP2(N(n), M(n), iters_p, iters_w, [], 'distributed', 'b', 0)];
         rms2 = [rms2; particle_filter_GP2(N(n), M(n), [], [], iters_c, 'distributed', 'c', 0)];
%        rms3 = [rms3; particle_filter_GPGA(N(n), iters_ga_g, 'gossip', [], 0)];
%         rms4 = [rms4; particle_filter_GPGA(N(n), iters_ga_b, 'broadcast', 'MFX', 0)];
%         rms5 = [rms5; particle_filter_GPFA(N(n), 1, iters_fa_g, 'gossip', 0)];
        fprintf('%d finished!\n',i);
    end
    armse(n,1) = mean(rms1,'all');
    armse(n,2) = mean(rms2,'all');
    armse(n,3) = mean(rms3,'all');
    armse(n,4) = mean(rms4,'all');
    armse(n,5) = mean(rms5,'all');
end
%save armse armse;
% 
%load armse;
% figure;
% plot(N, armse(:,1),'-s','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
% plot(N, armse(:,2),'-s','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
% plot(N, armse(:,3),'-s','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold on;
% plot(N, armse(:,4),'-s','Linewidth',1,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','MarkerSize',4); hold on;
% plot(N, armse(:,5),'-s','Linewidth',1,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerSize',4); hold off;
% figure;
% semilogy(N, armse(:,1),'-s','Linewidth',1,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize',4 ); hold on;
% semilogy(N, armse(:,2),'--s','Linewidth',1,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize',4); hold on;
% semilogy(N, armse(:,4),'-.^','Linewidth',1,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','MarkerSize',4); hold on;
% semilogy(N, armse(:,5),':h','Linewidth',1,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','MarkerSize',4); hold on;
% semilogy(N, armse(:,3),'-^','Linewidth',1,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','MarkerSize',4); hold off;
% legend('GPDPF(consensus)','GPDPF(diffusion)','GADPF(broadcast)','FADPF(gossip)','GADPF(gossip)')
% xlabel('particles'); ylabel('ARMSE');
%%
clear; clc; close all;
N_sensor = 9;
N = 50; % particle size
network = 'distributed';
%network = 'centralised';
iters = 10;
rmse = particle_filter_GP(N, network, iters, 1);




