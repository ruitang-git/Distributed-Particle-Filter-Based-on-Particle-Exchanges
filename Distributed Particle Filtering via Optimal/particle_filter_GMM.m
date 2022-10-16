function rmse = particle_filter_GMM(plt)

% ---------------------- %
% close all; clc; clear;
% plt = 0;
% configfile;
% trials = 20;
% 
% rmse_total = [];
% for trial = 1:trials
% ---------------------- %

% function rmse = particle_filter_GMM(plt,)
%
% Description:
% 'This is core function of the reproduction work.' 
%  this file reproduces the work of paper 
%
%  "Li, Jichuan, and Arye Nehorai. 
%  Distributed particle filtering via optimal fusion of Gaussian mixtures.
%  IEEE Transactions on Signal and Information Processing over Networks 4.2
%  (2017): 280-292."
%
% INPUTS:
%  plt - decide the plotting process
% 
% OUTPUTS:
% rmse - tracking error
%
configfile;
% [A,cd] = sensor(N_sensor,roi); % sensor location
[A, cd] = sensor_map(N_sensor, 200, 200, plt);
true_state = trajectory(STI,t,sigma_u);

% plot
if plt == 1
    figure;
    plot(true_state(1,:),true_state(2,:),'--o');
    hold on;
    gplot(A, cd.',':o');
    hold on;
    plot(cd(1,:),cd(2,:),'o');
    hold off;
    grid on;
    title('Agent trajectory');
end
    
x_ori = true_state(:,1); % origin state
d = length(x_ori); % state vector dimension
% % single sensor
% single_senor()  % function

x_est_out = []; % filtering output(estimated states)
rmse = [];

% Multiple sensors
A = A + eye(N_sensor); % make diagnol elements of ajacent matirx to 1
% Metroplis weight
D = sum(A,2); % degree(including itself)
EPS = zeros(N_sensor); % initialize
for i = 1:N_sensor
    for j = 1:N_sensor
        EPS(i,j) = 1/max(D(i),D(j));
    end
end
EPS = EPS .* (A - eye(N_sensor));
for k = 1:N_sensor
    EPS(k,k) = 1-sum(EPS(:,k));
end

%%
% filtering begin
% Global statics repository
DATA = {};
for sti = 1 : STI
    x = true_state(:,sti+1);
    
    % filtering
    for n_sensor = 1 : N_sensor
        l = cd(:,n_sensor); % certain sensor location
        if sti == 1 % particle initialization
            x_P = zeros(d,N); % particle set
            R = sigma_i^2 * eye(d);
            for n = 1 : N
                x_P(:,n) = x_ori + chol(R).' * randn(d,1);
            end
        else % generate particles 
            x_P = GMsample(DATA{n_sensor,4},DATA{n_sensor,5},DATA{n_sensor,6},N);
        end
        
        % observation with noise
        [h,R] = observation(x,l,sigma_v,sigma_w,1);
        % prior
        x_P_update = zeros(d,N);
        h_update = zeros(2,N);
        for n = 1 : N
            x_P_update(:,n) = transition(x_P(:,n),t,sigma_u);
            % observation without noise
            [h_update(:,n),~] = observation(x_P_update(:,n),l,sigma_v,sigma_w,0);
            P_w(n) = mvnpdf(h_update(:,n)', h', R);
        end
        % normalization.
        if sum(P_w) == 0
            P_w = 1/N*ones(1,N);
        else
            P_w = P_w./sum(P_w);
        end

        % Resampling
        % Me = 1/sum(P_w.^2); % efficient sample size: decide resampling
        Me = 0; % we first assume must resample which is easier for implementing
                % else we need to save the weights in previous time step
        if Me < 0.6*N
            for n = 1 : N
                x_P(:,n) = x_P_update(:,find(rand <= cumsum(P_w),1));   % particles with higher weights get more offsprings
            end
            P_w = ones(1,N)/N;
        else
            x_P = x_P_update;
        end
        
        % learn the prior prediction and proceed GMM approximation
        % save particles x_P_update
        GM_pk = fitgmdist(x_P_update',C,'RegularizationValue',1e-5,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
%         GM_pk = fitgmdist(x_P_update',C,'RegularizationValue',0.01,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
        a_pk = GM_pk.ComponentProportion';
        mu_pk = GM_pk.mu';
        Sigma_pk = diag_sigma(GM_pk.Sigma);
        Sigma_pk = reshape(Sigma_pk,[d,d*C]);
        % [a_pk, mu_pk, Sigma_pk] = GM(x_P_update,ones(N,1)/N,C); % have to decide how many components needed to be chosen ?
        % learn the local posterior and proceed GMM approximation
        % save particles x_P
        GM_k = fitgmdist(x_P',C,'RegularizationValue',1e-5,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
%         GM_k = fitgmdist(x_P',C,'RegularizationValue',0.01,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
        a_k = GM_k.ComponentProportion';
        mu_k = GM_k.mu';
        Sigma_k = diag_sigma(GM_k.Sigma);
        Sigma_k = reshape(Sigma_k,[d,d*C]);
        % [a_k, mu_k, Sigma_k] = GM(x_P,P_w',C); % have to decide how many components needed to be chosen ?
        
        DATA{n_sensor,1} = a_pk; DATA{n_sensor,2} = mu_pk; DATA{n_sensor,3} = Sigma_pk;
        DATA{n_sensor,4} = a_k; DATA{n_sensor,5} = mu_k; DATA{n_sensor,6} = Sigma_k;
    end
%  
%     % ---- plot(testing, delete in final version) @final sensor ---- %
%     figure;
%     plot(x(1),x(2),'kx','LineWidth',2,'MarkerSize',10);
%     hold on;
%     scatter(x_P_update(1,:),x_P_update(2,:),'b.');
%     hold on;
%     scatter(x_P(1,:),x_P(2,:),'r.');
%     hold off;
%     legend('true state', 'prior', 'posterior');
%     title('before GMM approximation');
%     % ---- plot(testing, delete in final version) @final sensor ---- %
    
%     % ---- verify the GMM approximation (testing, delete in final version) @final sensor ---- %
%     n_sensor =3;
%     p_pk = GMsample(DATA{n_sensor,1}, DATA{n_sensor,2}, DATA{n_sensor,3}, N);
%     p_k = GMsample(DATA{n_sensor,4}, DATA{n_sensor,5}, DATA{n_sensor,6}, N);
%     figure;
%     plot(x(3),x(4),'kx','LineWidth',2,'MarkerSize',10);
%     hold on;
%     scatter(p_pk(3,:),p_pk(4,:),'b.');
%     hold on;
%     scatter(p_k(3,:),p_k(4,:),'r.');
%     hold off;
%     legend('true state', 'prior', 'posterior');
%     title('after GMM approximation');
%     % ---- verify the GMM approximation(testing, delete in final version) @final sensor ---- %
    
    % ---- centralised -------%
%     centralised;
    % ---- centralised -------%
    
    % plot gaussian
    if plt == 1
        X=linspace(round(min(x_P(1,:))-2),round(max(x_P(1,:))+2),100);
        plot_GMM;
    end
    
    % fusion
    MEMORY = {}; % {{ },{ },{ }}
    MEMORY_prior = {};
    % convergence analysis
    cvg_errors = [];
    %[cvg_error(iter1) cvg_error(iter2) ... ]
    %[e(@s1,iter1) e(@s1,iter2) ... ]
    %[e(@s2,iter1) e(@s2,iter2) ... ]
    %[ ...         ...          ... ]
    % consensus
    ite = 1;
    consensus = consensus_ite;
    while consensus
        % communication
        for n_sensor = 1:N_sensor
            Nk = find(A(n_sensor,:));
            memory = {};
            for nk = 1:length(Nk)
                memory{nk,1} = DATA{Nk(nk),4};
                memory{nk,2} = DATA{Nk(nk),5};
                memory{nk,3} = DATA{Nk(nk),6};
            end
            MEMORY{1,n_sensor} = memory;
        end
        % update
        for n_sensor = 1:N_sensor
            [a, mu, Sigma] = GMFusion(MEMORY{1,n_sensor}, EPS(:,n_sensor), C, N_fusion, maxiter);
            % update
            DATA{n_sensor,4} = a;
            DATA{n_sensor,5} = mu;
            DATA{n_sensor,6} = Sigma;
        end
%         fprintf('[%d-th state interval] [%d-th consensus iteration] finished!\n', sti, ite);
        
        ite = ite + 1;
        consensus = consensus - 1;
%         cvg_analysis; cvg_errors = [cvg_errors cvg_error];
   
        % plot gaussian
        if plt == 1
            % plot_GMM;
        end
        
    end
    
    % fine tuning
    ite = 1;
    f_tuning = f_tuning_ite;
    while f_tuning
        % communication
        for n_sensor = 1:N_sensor
            Nk = find(A(n_sensor,:));
            memory = {};
            % prior
            memory_prior = {};
            for nk = 1:length(Nk)
                % posterior
                memory{nk,1} = DATA{Nk(nk),4};
                memory{nk,2} = DATA{Nk(nk),5};
                memory{nk,3} = DATA{Nk(nk),6};
                % prior
                memory_prior{nk,1} = DATA{Nk(nk),1};
                memory_prior{nk,2} = DATA{Nk(nk),2};
                memory_prior{nk,3} = DATA{Nk(nk),3};
            end
            % switching
            % make the data of center sensor shifting to the first place
            pos = find(Nk == n_sensor);
            % posterior
            tmp_a = memory{pos,1};
            tmp_mu = memory{pos,2};
            tmp__sigma = memory{pos,3};
            memory{pos,1} = memory{1,1};
            memory{pos,2} = memory{1,2};
            memory{pos,3} = memory{1,3};
            memory{1,1} = tmp_a;
            memory{1,2} = tmp_mu;
            memory{1,3} = tmp__sigma;
            MEMORY{1,n_sensor} = memory;
            % prior
            tmp_a = memory_prior{pos,1};
            tmp_mu = memory_prior{pos,2};
            tmp__sigma = memory_prior{pos,3};
            memory_prior{pos,1} = memory_prior{1,1};
            memory_prior{pos,2} = memory_prior{1,2};
            memory_prior{pos,3} = memory_prior{1,3};
            memory_prior{1,1} = tmp_a;
            memory_prior{1,2} = tmp_mu;
            memory_prior{1,3} = tmp__sigma;
            MEMORY_prior{1,n_sensor} = memory_prior;
        end
        % update
        for n_sensor = 1:N_sensor
            % posterior
            [a, mu, Sigma] = fine_tuning(MEMORY{1,n_sensor});
            % update
            DATA{n_sensor,4} = a;
            DATA{n_sensor,5} = mu;
            DATA{n_sensor,6} = Sigma;
            % prior
            [a, mu, Sigma] = fine_tuning(MEMORY_prior{1,n_sensor});
            % update
            DATA{n_sensor,1} = a;
            DATA{n_sensor,2} = mu;
            DATA{n_sensor,3} = Sigma;
        end
%         fprintf('[%d-th state interval] [%d-th fine-tuning iteration] finished!\n', sti, ite);

        ite = ite + 1;
        f_tuning = f_tuning - 1;
%         cvg_analysis; cvg_errors = [cvg_errors cvg_error];
    end
    
    % plot gaussian
    if plt == 1
        plot_GMM;
    end

    % average convergence error
    mean_cvg_errors = mean(cvg_errors,1);
    %
    if plt == 1
        figure;
        for n_sensor = 1:N_sensor
            plot(cvg_errors(n_sensor,:),'-o');
            hold on;
        end
        plot(mean_cvg_errors,'--*');
        hold off;
        title('KL Distance');
        legend('sensor1','sensor2','sensor3','sensor4','sensor5','sensor6','sensor7','sensor8','sensor9','average');
    end

    % ---- verify the fusion(testing, delete in final version) @1_st sensor ---- %
%     p_1 = GMsample(DATA{1,4}, DATA{1,5}, DATA{1,6}, N);
%     p_2 = GMsample(DATA{1,1}, DATA{1,2}, DATA{1,3}, N);
%     figure;
%     plot(x(1),x(2),'kx','LineWidth',2,'MarkerSize',10);
%     hold on;
%     scatter(p_1(1,:),p_1(2,:),'r.');
%     hold on;
%     scatter(p_2(1,:),p_2(2,:),'b.');
%     hold off;
%     legend('true state', 'posterior product','prior');
%     title('GMM approximation fusion');
    % ---- verify the fusion(testing, delete in final version) @1_st sensor ---- %
        
    % recovery
    for n_sensor = 1:N_sensor
        data = DATA(n_sensor,:);
        % Global posterior
        % [a, Mu, Sigma] = GMRecovery(x, data, K, C, M, maxiter, reg, plt)
        [a, Mu, Sigma] = GMRecovery(x, data, N_sensor, C, N_recovery, maxiter,1e-5, 0);
        DATA{n_sensor,4} = a;
        DATA{n_sensor,5} = Mu;
        DATA{n_sensor,6} = Sigma;
%         plt = 0;
    end
    
    % ---- verify the recovery(testing, delete in final version) @1_st sensor ---- %
%     p_1 = GMsample(DATA{3,4}, DATA{3,5}, DATA{3,6}, N);
%     figure;
%     plot(x(1),x(2),'kx','LineWidth',2,'MarkerSize',10);
%     hold on;
%     scatter(p_1(1,:),p_1(2,:),'r.');
%     hold off;
%     legend('true state', 'global posterior');
%     title('after recovery');
    % ---- verify the recovery(testing, delete in final version) @1_st sensor ---- %
        
    % calculate estimated trajectory
    % use the data at the first sensor as a reference
    x_est = [];
    for n_sensor = 1:N_sensor
        x_est = [x_est DATA{n_sensor,5}*DATA{n_sensor,4}];
    end
    x_est_out = [x_est_out mean(x_est,2)];
    error = 0;
    for n_sensor = 1:N_sensor
        error = error + norm(x_est(:,n_sensor)-x);
    end
    rmse = [rmse error/N_sensor];
end

if plt == 1
    figure;
    p1 = plot(true_state(1,2:end),true_state(2,2:end),'--o');
    hold on;
    p2 = plot(x_est_out(1,:),x_est_out(2,:),'--or'); % excludes the intial state estimate
    hold on;
    gplot(A, cd.',':o');
    hold off;
    grid on;
    legend([p1,p2],{'true','estimated'});
    title('Agent trajectory');
%     figure;
%     plot(rmse);
%     title('RMSE Error');
%     xlabel('Time');ylabel('RMSE');
end

% -------------------- %
% rmse_total = [rmse_total; rmse];
% fprintf('%d iter finished!\n', trial);
end
% -------------------- %
