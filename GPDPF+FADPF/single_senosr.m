%% SPF
% function single_senosr()

%%% result
% n=50;350
% n=200;200+
% n=500;160+

clear; clc;
close all;

config;

% Monte Carlo
MC = 100;
rmse_wp_global = cell(MC, 1);
Me_global = cell(MC,1);
% parfor mc = 1:MC
for mc = 1:MC
tic;    
%     sigma_u = 0.1;
%     sigma_u = 0.5;
%     STI = 50; % state transition intervals
%     N = 50; % particle size

    rmse_wp = zeros(1,1+STI);
    
%     rng(0);
    true_state = trajectory(STI,t,sigma_u);
%     rng('shuffle');
    sl = [0;0]; % single sensor
    x = true_state(:,1); % origin state vector
    d = length(x);
    x_P = zeros(d,N); % particle set
%     R = [0.1 0 0 0; 0 0.1 0 0; 0 0 0.01 0; 0 0 0 0.01];
    R = [0.1 0 0 0 0 0; 0 0.1 0 0 0 0; 0 0 0.01 0 0 0; 0 0 0 0.01 0 0; 0 0 0 0 0.01 0; 0 0 0 0 0 0.01];
    for n = 1 : N
        x_P(:,n) = x + chol(R).' * rand(d,1);
    end
    P_w = ones(1,N)/N;

    x_est = mean(x_P,2);
    x_est_out_wp = x_est; % the vector of particle filter estimates.
    rmse_wp(1) = rmse_wp(1) + norm(x_est(1:2) - x(1:2))^2; % RMS position error
    
    % particle efficient size
    Me = zeros(1,STI);
    for sti = 1 : STI
        x = true_state(:,sti+1);
        % observation with noise
        [h,R] = observation(x,sl,1);
        % prior
        x_P_update = zeros(d,N);
        h_update = zeros(2,N);
        for n = 1 : N
            % x_P_update(:,n) = transition_cv(x_P(:,n),sigma_u); % CV model is used
            x_P_update(:,n) = transition(x_P(:,n),t,sigma_u);
            % observation without noise
            [h_update(:,n),~] = observation(x_P_update(:,n),sl,0);
            P_w(n) = P_w(n) * mvnpdf(h_update(:,n)', h', R);
        end
        % normalization.
        P_w = max(P_w, 1e-10);
        P_w = P_w./sum(P_w);
        
        %%%%% state estimate
        % a) weigthed particles
        x_est_wp = x_P_update*P_w';

        %%%% resampling
%         Me = 0;
        Me(sti) = 1/sum(P_w.^2); % efficient sample size: decide resampling | must resample
%         Me(sti) = 0;
        if Me(sti) < 0.6*N
            % 1) standard resampling method
            [ x_P, P_w ] = Resample( x_P_update, P_w );
        else
            x_P = x_P_update;
        end

        error_wp = norm(x(1:2) - x_est_wp(1:2))^2;
        rmse_wp(sti+1) = rmse_wp(sti+1) + error_wp;
        x_est_out_wp = [x_est_out_wp x_est_wp]; % estimated state

    end

    rmse_wp_global{mc, 1} = rmse_wp; 
    Me_global{mc,1} = Me;

    if mod(mc,10) == 0
        figure;
        p1 = plot(true_state(1,:),true_state(2,:),'-b');
        hold on;
        p2 = plot(x_est_out_wp(1,:),x_est_out_wp(2,:),'-r','LineWidth',1);
        hold off;
        title('Agent trajectory');
    end
toc;
end

rmse_wp = sqrt(sum(cell2mat(rmse_wp_global),1)/(2*MC));
Me = sum(cell2mat(Me_global),1)/MC;


figure;
plot(rmse_wp,'--');
xlabel('time/s'); ylabel('rmse');
figure;
plot(log(rmse_wp),'--');
xlabel('time/s'); ylabel('log(rmse)');

figure;
plot(Me,'--');
xlabel('time/s'); ylabel('Effective sample size'); grid on;
%% single sensor
% function single_senosr()
clear; clc;
close all;

config;

% Monte Carlo
MC = 5;
p_var = [0 0.2 0.4 0.6];
rmse_wp_global = cell(MC, 1);
rmse_mmse_global = cell(MC, 1);
rmse_map_global = cell(MC, 1);
Me_global = cell(MC,1);
% parfor mc = 1:MC
for mc = 1:MC
tic;    
%     sigma_u = 0.1;
%     sigma_u = 0.5;
%     STI = 50; % state transition intervals
%     N = 50; % particle size

    rmse_wp = zeros(1,1+STI);
    rmse_mmse = zeros(1,1+STI);
    rmse_map = zeros(1,1+STI);
    
%     rng(0);
    true_state = trajectory(STI,t,sigma_u);
%     rng('shuffle');
    sl = [0;0]; % single sensor
    x = true_state(:,1); % origin state vector
    d = length(x);
    x_P = zeros(d,N); % particle set
%     R = [0.1 0 0 0; 0 0.1 0 0; 0 0 0.01 0; 0 0 0 0.01];
    R = [0.1 0 0 0 0 0; 0 0.1 0 0 0 0; 0 0 0.01 0 0 0; 0 0 0 0.01 0 0; 0 0 0 0 0.01 0; 0 0 0 0 0 0.01];
    for n = 1 : N
        x_P(:,n) = x + chol(R).' * rand(d,1);
    end
    P_w = ones(1,N)/N;

    x_est = mean(x_P,2);
    x_est_out_wp = x_est; % the vector of particle filter estimates.
    x_est_out_mmse = x_est;
    x_est_out_map = x_est;
    rmse_wp(1) = rmse_wp(1) + norm(x_est(1:2) - x(1:2))^2; % RMS position error
    rmse_mmse(1) = rmse_mmse(1) + norm(x_est(1:2) - x(1:2))^2;
    rmse_map(1) = rmse_map(1) + norm(x_est(1:2) - x(1:2))^2;
    
    % particle efficient size
    Me = zeros(1,STI);
    for sti = 1 : STI
        x = true_state(:,sti+1);
        % observation with noise
        [h,R] = observation(x,sl,1);
        % prior
        x_P_update = zeros(d,N);
        h_update = zeros(2,N);
        for n = 1 : N
            % x_P_update(:,n) = transition_cv(x_P(:,n),sigma_u); % CV model is used
            x_P_update(:,n) = transition(x_P(:,n),t,sigma_u);
            % observation without noise
            [h_update(:,n),~] = observation(x_P_update(:,n),sl,0);
            P_w(n) = P_w(n) * mvnpdf(h_update(:,n)', h', R);
        end
        % normalization.
        P_w = max(P_w, 1e-10);
        P_w = P_w./sum(P_w);
        % particle interpolation
        % generate a new partlce set with size N*(N-1)
        % this is the preparation step for GP-based resampling
        
        
        particle = x_P_update.'; % old set
        particle_set = particle;
        N = size(particle,1); % particle: N*d
        for i = 1:N-1
            new_particle = particle(i,:) + (particle(i+1,:)-particle(i,:))/2;
            particle_set = [particle_set; new_particle];
        end


        %%%%% GP 
        %%% Do I have to choose the lower and upper bound for
        %%% hyperparamters ???
%         tic;
%         [alpha, l, sigma_n, y_hat] = GP(particle, P_w',x, particle_set); % y_hat: y+3*sigma; required in the resamping stage
        [alpha, l, sigma_n, y_hat] = GP_optfree(particle, P_w',x, particle_set); % y_hat: y+3*sigma; required in the resamping stage
        y_hat = y_hat/sum(y_hat);
%         toc;

        %%%%% state estimate
        % a) weigthed particles
        x_est_wp = x_P_update*P_w';
        % b) MMSE
        x_est_mmse = x_P_update*alpha/sum(alpha);
        % c) MAP
        [~,index] = max(P_w);
        inits = particle(index,:);
        x_est_map = fmincon(@(x_est_map) MAPobj(x_est_map,particle,alpha,l,sigma_n), inits,[], [], [], [], min(particle,[],1), max(particle,[],1));
        x_est_map = x_est_map.';

        %%%% resampling
%         Me = 0;
        Me(sti) = 1/sum(P_w.^2); % efficient sample size: decide resampling | must resample
%         Me(sti) = 0;
        if Me(sti) < 0.6*N
            % 1) standard resampling method
            [ x_P, P_w ] = Resample( x_P_update, P_w );
%             debug = 1;
            [ x_P_res, P_w_res ] = Resample( x_P_update, P_w );
            % 2) GP-based resampling method
            % step1) GP approx.
%             tic;
%             [alphaf, lf, sigma_nf, ~] = GP(particle_set, y_hat, x);
            [alphaf, lf, sigma_nf, ~] = GP_optfree(particle_set, y_hat, x);
%             toc;
            % step2) create GMM
            sigma = l^2 * ones(1,d); p = sqrt(2*pi)*sigma_n*alpha/((sigma_n*sqrt(2*pi))^d*sum(alpha)); p = max(p, 1e-5); p = p/sum(p);
            sigmaf = lf^2 * ones(1,d); pk = sqrt(2*pi)*sigma_nf*alphaf/((sigma_nf*sqrt(2*pi))^d*sum(alphaf)); pk = max(pk, 1e-5); pk = pk/sum(pk);
            % step3) sampling 
            for n = 1:N
                u =  rand; 
                if u < p_var(1)
                    gm = gmdistribution(particle_set, sigmaf, pk); % based on alphaf
                else
                    gm = gmdistribution(particle, sigma, p); % based on alpha
                end
                x_P(:,n) = random(gm,1).';
            end
            % plot
            % figure;plot(x(1),x(2),'kx','LineWidth',2,'MarkerSize',10); hold on;plot(x_P(1,:),x_P(2,:),'b.');hold on; plot(x_P_update(1,:),x_P_update(2,:),'r.'); hold off;
            P_w = ones(1,N)/N;
        else
            x_P = x_P_update;
        end

        
%         gprMdl = fitrgp(x_P(:,:)', P_w, 'Basis', 'linear', 'FitMethod','exact', 'PredictMethod','exact');
%         figure;plot3(x_P(1,:),x_P(2,:), P_w ,'b.'); hold on; plot3(x_P(1,:),x_P(2,:),ypred,'r.'); hold off;
%         x_cord = min(x_P(1,:)):0.1:max(x_P(1,:)); y_cord = min(x_P(2,:)):0.1:max(x_P(2,:));
%         [X,Y] = meshgrid(x_cord,y_cord);
%         cord = [];
%         for i = 1:size(X,1)
%             for j = 1:size(X,2)
%                 cord = [cord; [X(i,j) Y(i,j)]];
%             end
%         end
%         ypred = predict(gprMdl,cord);
%         figure;plot3(x_P(1,:),x_P(2,:), P_w ,'b.'); hold on; surf(X,Y,(reshape(ypred,[size(X,2), size(X,1)]))'); hold off;
%         SigmaL = gprMdl.KernelInformation.KernelParameters(1);
%         SigmaF = gprMdl.KernelInformation.KernelParameters(2);
%         K = zeros(size(x_P,2)); 
%         for i = 1:size(K,1)
%             for j = 1:size(K,2)
%                 x1 = x_P(:,i); x2 = x_P(:,j);
%                 K(i,j) = SigmaF^2 * exp(-(x1-x2)'*(x1-x2)/(2*SigmaL^2));
%             end
%         end
%         alpha = max(inv(K) * P_w',0);

        error_wp = norm(x(1:2) - x_est_wp(1:2))^2;
        rmse_wp(sti+1) = rmse_wp(sti+1) + error_wp;
        x_est_out_wp = [x_est_out_wp x_est_wp]; % estimated state
        error_mmse = norm(x(1:2) - x_est_mmse(1:2))^2;
        error_map = norm(x(1:2) - x_est_map(1:2))^2;
        rmse_mmse(sti+1) = rmse_mmse(sti+1) + error_mmse;
        rmse_map(sti+1) = rmse_map(sti+1) + error_map;
        x_est_out_mmse = [x_est_out_mmse x_est_mmse];
        x_est_out_map = [x_est_out_map x_est_map];
    end

    rmse_wp_global{mc, 1} = rmse_wp;
    rmse_mmse_global{mc, 1} = rmse_mmse; 
    rmse_map_global{mc, 1} = rmse_map;   
    Me_global{mc,1} = Me;

    if mod(mc,1) == 0
        figure;
        p1 = plot(true_state(1,:),true_state(2,:),'-b');
        hold on;
        p2 = plot(x_est_out_wp(1,:),x_est_out_wp(2,:),'-r','LineWidth',1);
        hold on;
        p3 = plot(x_est_out_mmse(1,:),x_est_out_mmse(2,:),'-g','LineWidth',1);
        hold on;
        p4 = plot(x_est_out_map(1,:),x_est_out_map(2,:),'-y','LineWidth',1);
        hold off;
        grid on;
        legend([p1,p2,p3,p4],{'true','estimated(weighted particles)','GP-MMSE','GP-MAP'});
        title('Agent trajectory');
    end
toc;
end

rmse_wp = sqrt(sum(cell2mat(rmse_wp_global),1)/(2*MC));
rmse_mmse = sqrt(sum(cell2mat(rmse_mmse_global),1)/(2*MC));
rmse_map = sqrt(sum(cell2mat(rmse_map_global),1)/(2*MC));
Me = sum(cell2mat(Me_global),1)/MC;

save rmse_wp_global rmse_wp_global;
save rmse_mmse_global  rmse_mmse_global;
save rmse_map_global rmse_map_global;
save rmse_wp rmse_wp;
save rmse_mmse rmse_mmse;
save rmse_map rmse_map;

figure;
plot(rmse_wp,'--');
hold on;
plot(rmse_mmse,'--');
hold on;
plot(rmse_map,'--');
hold off;
legend('SPF','GP-MMSE','GP-MAP');
xlabel('time/s'); ylabel('rmse');
figure;
plot(log(rmse_wp),'--');
hold on;
plot(log(rmse_mmse),'--');
hold on;
plot(log(rmse_map),'--');
hold off;
legend('SPF','GP-MMSE','GP-MAP');
xlabel('time/s'); ylabel('log(rmse)');

figure;
plot(Me,'--');
xlabel('time/s'); ylabel('Effective sample size'); grid on;

%%
function [ particle_res, weight_res ] = Resample( particle, weight )
%RESAMPLE 이 함수의 요약 설명 위치
%   자세한 설명 위치
[r, numParticle] = size(particle);

particle_tmp = zeros(r, numParticle);

weight = weight/sum(weight);
weight_cdf = cumsum(weight);
for j = 1:1:numParticle
    index_find = find(rand <= weight_cdf, 1);
    if isempty(index_find)
        [a, index_find] = max(weight);
    end
    particle_tmp(:,j) = particle(:,index_find);
end
particle_res = particle_tmp;
weight_res = ones(1,numParticle)/numParticle;

end



function fcn = MAPobj (params,X,alpha,l,sigman)
n = length(alpha);
for i = 1:n
    kstar(i) = exp(-(norm(X(i,:)-params))^2/(2*l^2));
    if X(i,:) == params
        kstar(i) = kstar(i) + sigman^2;
    end
end
fcn = -kstar*alpha; 
end


