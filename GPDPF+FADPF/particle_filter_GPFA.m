function rmse = particle_filter_GPFA(N, iters_FA, iters_com, protocol, plt)

% Input
% N: particle size
% iters_FA: 1(no use)
% iters_com: number of communication iterations
% protocol: 'gossip'
% plt: plot figure or not 

% clear all; close all; clc;
% protocol = 'broadcast';
% % protocol = 'gossip';
% N = 5; 
% iters_GA = 10;
%

t=1;
sigma_u = 0.5;
sigma_v = 1;
sigma_w = 1;
sigma_i = 1; % particle initialization
STI = 30; % state transition intervals
N_sensor = 9; % number of sensors
[A, cd] = sensor_map(N_sensor, 200, 200, plt);
true_state = trajectory(STI,t,sigma_u);
p_var = [0 0.2 0.4 0.6];
    
x_ori = true_state(:,1); % origin state
d = length(x_ori); % state vector dimension

x_est_out = cell(N_sensor,1); % filtering output(estimated states)
rmse = [];

%%
x_P = cell(N_sensor,1); P_w = cell(N_sensor,1); particle = cell(N_sensor,1);

for n_sensor = 1:N_sensor
    x_P{n_sensor,1} = zeros(d,N); % particle set
    P_w{n_sensor,1} = ones(1,N)/N;
    R = sigma_i^2 * eye(d);
    for n = 1 : N
        x_P{n_sensor,1}(:,n) = x_ori + chol(R).' * randn(d,1);
    end
end

for sti = 1 : STI
    x = true_state(:,sti+1); 
    h_global = cell(N_sensor,1);
    neighbors = {[2 4], [1 3 5], [2 6], [1 5 7], [2 4 6 8], [3 5 9], [4 8], [5 7 9], [6 8]}; % 9-node grid network 
    for n_sensor = 1 : N_sensor
        ls = cd(:,n_sensor); % certain sensor location

        % observation with noise
        [h,R] = observation(x,ls,1);
        h_global{n_sensor,1} = h;
        % prior
        % x_P_update = zeros(d,N);
        h_update = zeros(2,N);

        for n = 1 : N
            x_P{n_sensor,1}(:,n) = transition(x_P{n_sensor,1}(:,n),t,sigma_u);
            % observation without noise
            [h_update(:,n),~] = observation(x_P{n_sensor,1}(:,n),ls,0);
            P_w{n_sensor,1}(n) = P_w{n_sensor,1}(n) * mvnpdf(h_update(:,n)', h', R);
        end

        % normalization.
        P_w{n_sensor,1} = max(P_w{n_sensor,1}, 1e-10);
        P_w{n_sensor,1} = P_w{n_sensor,1}./sum(P_w{n_sensor,1});

        %%%%% GP-enhanced density estimation and resampling
        particle{n_sensor,1} = x_P{n_sensor,1}.'; % old set
        particle_set = particle{n_sensor,1};
%             for i = 1:N-1
%                 new_particle = particle{n_sensor,1}(i,:) + (particle{n_sensor,1}(i+1,:)-particle{n_sensor,1}(i,:))/2;
%                 particle_set = [particle_set; new_particle];
%             end
        [alpha, l, sigma_n, ~] = GP_optfree(particle{n_sensor,1}, P_w{n_sensor,1}',x);
%             [alpha, l, sigma_n, y_hat] = GP_optfree(particle{n_sensor,1}, P_w{n_sensor,1}',x, particle_set); % y_hat: y+3*sigma; required in the resamping stage
%             y_hat = y_hat/sum(y_hat);
%             [alphaf, lf, sigma_nf, ~] = GP_optfree(particle_set, y_hat, x);
        sigma = l^2 * ones(1,d); p = sqrt(2*pi)*sigma_n*alpha/((sigma_n*sqrt(2*pi))^d*sum(alpha)); p = max(p, 1e-5); p = p/sum(p);
%             sigmaf = lf^2 * ones(1,d); pk = sqrt(2*pi)*sigma_nf*alphaf/((sigma_nf*sqrt(2*pi))^d*sum(alphaf)); pk = max(pk, 1e-5); pk = pk/sum(pk);
        x_P{n_sensor,1} = [];
        for n = 1:N
            u =  rand; 
            if u < p_var(1)
                gm = gmdistribution(particle_set, sigmaf, pk); % based on alphaf
            else
                gm = gmdistribution(particle{n_sensor,1}, sigma, p); % based on alpha
            end
            x_P{n_sensor,1}(:,n) = random(gm,1).';
        end 
        
        %%% standard resampling
%         [ x_P{n_sensor,1}, ~ ] = Resample( x_P{n_sensor,1}, P_w{n_sensor,1}, N );

    end
    %%%%% fA fusion
    x_P = FA(x_P, cd, h_global, R, iters_FA, iters_com, neighbors, protocol, x);
    for n_sensor = 1:N_sensor
         P_w{n_sensor,1} = ones(1,N)/N;
    end
   
    % calculate estimated trajectory
    % use the data at the first sensor as a reference
    for n_sensor = 1:N_sensor
        x_est = mean(x_P{n_sensor,1},2); x_est_out{n_sensor,1} = [x_est_out{n_sensor,1} x_est];
    end
    error = 0;
    for n_sensor = 1:N_sensor
        x_est = mean(x_P{n_sensor,1},2);
        error = error + norm(x_est-x);
    end
    rmse = [rmse error/N_sensor];
end

if plt == 1
    for n_sensor = 1:1
        figure;
        p1 = plot(true_state(1,2:end),true_state(2,2:end),'--o'); hold on;
        p2 = plot(x_est_out{n_sensor,1}(1,:),x_est_out{n_sensor,1}(2,:),'--or'); hold on;% excludes the intial state estimate
        gplot(A, cd.',':o'); hold off; grid on;
        legend([p1,p2],{'true','estimated'}); title('Agent trajectory');
    end
end
%     figure;
%     plot(rmse);
%     title('RMSE Error');
%     xlabel('Time');ylabel('RMSE');

function [ particle_res, weight_res ] = Resample( particle, weight, numParticle )
%RESAMPLE 이 함수의 요약 설명 위치
%
if nargin < 3
    [r, numParticle] = size(particle);
else
    r = size(particle,1);
end

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

end
