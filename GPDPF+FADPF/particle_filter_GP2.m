function rmse = particle_filter_GP2(N, M, iters_p, iters_w, iters, network, distributed, plt)
% input
% N: local particle size
% M: particle size in unit center  = round(N*phi*N_sensor)
% iters_p: consensus iterations on particles; (work if distributed = 'b') 
% iters_w: consensus iterations on weights; (work if distributed = 'b') 
% iters: diffusion iterations (work if distributed = 'c') 
% network: 'distributed'
% distributed: 'b'(consensus-based) OR 'c'(diffusion-based)
% plt: plot or not

t=1;
sigma_u = 0.5;
sigma_v = 1;
sigma_w = 1;
sigma_i = 1; % particle initialization
STI = 30; % state transition intervals
N_sensor = 9; % number of sensors
roi = 120; % parameter determining the connectivity
Sigma = sigma_u^2 * ...
     [1/20*t^5 0 1/8*t^4 0 1/6*t^3 0;...
      0 1/20*t^5 0 1/8*t^4 0 1/6*t^3;...
      1/8*t^4 0 1/3*t^3 0 1/2*t^2  0;...
      0 1/8*t^4 0 1/3*t^3 0 1/2*t^2 ;...
      1/6*t^3 0 1/2*t^2 0 t 0       ;...
      0 1/6*t^3 0 1/2*t^2 0 t      ];

% [A,cd] = sensor(N_sensor,roi); % sensor location
[A, cd] = sensor_map(N_sensor, 200, 200, plt);
true_state = trajectory(STI,t,sigma_u);
p_var = [0 0.2 0.4 0.6];


% plot
% if plt == 1
%     figure;
%     plot(true_state(1,:),true_state(2,:),'--o'); hold on;
%     gplot(A, cd.',':o'); hold on;
%     plot(cd(1,:),cd(2,:),'o'); hold off;
%     grid on; title('Agent trajectory');
% end
    
x_ori = true_state(:,1); % origin state
d = length(x_ori); % state vector dimension
% % single sensor
% single_senor()  % function

x_est_out = cell(N_sensor,1); % filtering output(estimated states)
rmse = [];

%%
% filtering begin
% Global statics repository
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
    
    %%%%% filtering
    if strcmp(network,'centralised')
        % 1) fusion center based
        h_global = cell(N_sensor,1);
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
%             particle_set = particle{n_sensor,1};
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
            for n = 1:N
                u =  rand; 
                if u < p_var(1)
                    gm = gmdistribution(particle_set, sigmaf, pk); % based on alphaf
                else
                    gm = gmdistribution(particle{n_sensor,1}, sigma, p); % based on alpha
                end
                x_P{n_sensor,1}(:,n) = random(gm,1).';
            end 
%             [ x_P{n_sensor,1}, ~ ] = Resample( x_P{n_sensor,1}, P_w{n_sensor,1} );
     
        end
        
        % communication (4 steps)
        % step1) PE - CU (transmitted: local particles)
        x_P_center = x_P; P_w_center = cell(N_sensor,1);
        for n_sensor = 1:N_sensor
            x_P_center{n_sensor,1} = x_P_center{n_sensor,1}(:,randperm(N, floor(M/N_sensor)));
        end
        % step2) CU - PE (transmitted: particles)
        % step3) PE - CU (transmitted: weights)
        for n_sensor = 1:N_sensor
            xstar = x_P_center{n_sensor,1}'; % 50x6
            for nn_sensor = 1:N_sensor
                if nn_sensor ~= n_sensor
                    % 
                    bestEstimate = zeros(size(xstar,1),1);
                    for n = 1:size(xstar,1)
                        [hp,~] = observation(xstar(n,:)',cd(:,nn_sensor),0);
                        bestEstimate(n) = mvnpdf(hp', h_global{nn_sensor,1}', R);
                    end
                    % (GP)
%                     X = particle{nn_sensor,1};
%                     y = P_w{nn_sensor,1}';
%                     [~, ~, ~, bestEstimate, ~] = GP_optfree(X, y, x, xstar);
                    P_w_center{n_sensor,1} = [P_w_center{n_sensor,1} bestEstimate];
                end
            end
            % sum
            P_w_center{n_sensor,1} = sum(P_w_center{n_sensor,1}, 2); 
            % product
%             tmp = ones(size(P_w_center{n_sensor,1},1),1);
%             for n = 1:size(P_w_center{n_sensor,1},2)
%                 tmp = tmp.*P_w_center{n_sensor,1}(:,n);
%             end
%             P_w_center{n_sensor,1} = tmp;            
        end
        P_w_center = cell2mat(P_w_center); P_w_center = reshape(P_w_center,[],1);
        x_P_center = cell2mat(x_P_center');
        [ x_P_center, ~ ] = Resample( x_P_center, P_w_center, N );
        % step4) CU - PE (transmitted: particles)
        for n_sensor = 1:N_sensor % backward (PE -> Center)
            x_P{n_sensor,1} = x_P_center;
            P_w{n_sensor,1} = ones(1,N)/N;
        end
  
    elseif strcmp(network,'distributed')
        if strcmp(distributed, 'a') % GP is used to indicate the fitness
            % 2.a) distributed resampling
            neighbors = {[2 4], [1 3 5], [2 6], [1 5 7], [2 4 6 8], [3 5 9], [4 8], [5 7 9], [6 8]}; % 9-node grid network 
            for n_sensor = 1 : N_sensor
                ls = cd(:,n_sensor); % certain sensor location

                % observation with noise
                [h,R] = observation(x,ls,1);
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
                for i = 1:N-1
                    new_particle = particle{n_sensor,1}(i,:) + (particle{n_sensor,1}(i+1,:)-particle{n_sensor,1}(i,:))/2;
                    particle_set = [particle_set; new_particle];
                end
                [alpha, l, sigma_n, y_hat] = GP_optfree(particle{n_sensor,1}, P_w{n_sensor,1}',x, particle_set); % y_hat: y+3*sigma; required in the resamping stage
                y_hat = y_hat/sum(y_hat);
                [alphaf, lf, sigma_nf, ~] = GP_optfree(particle_set, y_hat, x);
                sigma = l^2 * ones(1,d); p = sqrt(2*pi)*sigma_n*alpha/((sigma_n*sqrt(2*pi))^d*sum(alpha)); p = max(p, 1e-5); p = p/sum(p);
                sigmaf = lf^2 * ones(1,d); pk = sqrt(2*pi)*sigma_nf*alphaf/((sigma_nf*sqrt(2*pi))^d*sum(alphaf)); pk = max(pk, 1e-5); pk = pk/sum(pk);
                for n = 1:N
                    u =  rand; 
                    if u < p_var(1)
                        gm = gmdistribution(particle_set, sigmaf, pk); % based on alphaf
                    else
                        gm = gmdistribution(particle{n_sensor,1}, sigma, p); % based on alpha
                    end
                    x_P{n_sensor,1}(:,n) = random(gm,1).';
                end 
            end
            %%%%% a) communication (gossip is used)
            iters = 20;
            for iter = 1:iters
                id1 = randperm(N_sensor,1); id2 = neighbors{id1}(randperm(length(neighbors{id1}),1));  
                xstar = x_P{id2,1}';
                X = particle{id1,1}; y = P_w{id1,1}';
                [~, ~, ~, bestEstimate, ~] = GP_optfree(X, y, x, xstar);
                P_w{id1,1} = bestEstimate';
                xstar = x_P{id1,1}';
                X = particle{id2,1}; y = P_w{id2,1}';
                [~, ~, ~, bestEstimate, ~] = GP_optfree(X, y, x, xstar);
                P_w{id2,1} = bestEstimate';
                P_w{id1,1} = [P_w{id2,1} P_w{id1,1}]; 
                x_P{id1,1} = [x_P{id1,1} x_P{id2,1}];
                [ x_P{id1,1}, P_w{id1,1} ] = Resample( x_P{id1,1}, P_w{id1,1}, N );
                [ x_P{id2,1}, P_w{id2,1} ] = Resample( x_P{id1,1}, P_w{id1,1}, N );
            end  

            %%%% b) communication (group based (adjacent nodes))
            iters = 20;
            for iter = 1:iters
                id1 = randperm(N_sensor,1); id2 = neighbors{id1};
                for id = 1:length(id2)
                    xstar = x_P{id1,1}';
                    X = particle{id2(id),1}; y = P_w{id2(id),1}';
                    [~, ~, ~, bestEstimate, ~] = GP_optfree(X, y, x, xstar);
                    P_w{id2(id),1} = [ones(1,N)/N bestEstimate'];
                    P_w{id1,1} = [P_w{id1,1}; bestEstimate'];
                    x_P{id2(id),1} = [x_P{id2(id),1} x_P{id1,1}];
                    [ x_P{id2(id),1}, P_w{id2(id),1} ] = Resample( x_P{id2(id),1}, P_w{id2(id),1}, N );
                end
                [ x_P{id1,1}, P_w{id1,1} ] = Resample( x_P{id1,1}, sum(P_w{id1,1}), N );
            end  
            
        elseif strcmp(distributed, 'b')
            % 2.b) distributed resampling
            h_global = cell(N_sensor,1);
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
                for n = 1:M/N_sensor
                    u =  rand; 
                    if u < p_var(1)
                        gm = gmdistribution(particle_set, sigmaf, pk); % based on alphaf
                    else
                        gm = gmdistribution(particle{n_sensor,1}, sigma, p); % based on alpha
                    end
                    x_P{n_sensor,1}(:,n) = random(gm,1).';
                end 

            end
            %%%%% communication (particle exchange(consensus is used))
            data_particles = [];
            data_weights = [];
            for n_sensor = 1:N_sensor % communication can be reduced using sparse matrix
                tmp = x_P{n_sensor,1};
                x_P{n_sensor,1} = zeros(d,M);
                x_P{n_sensor,1}(:,(n_sensor-1)*M/N_sensor+1:n_sensor*M/N_sensor) = tmp;
                x_P{n_sensor,1}= reshape(x_P{n_sensor,1},[],1);
                data_particles = [data_particles x_P{n_sensor,1}];
            end
            ave_x_P = gossip(A, data_particles, iters_p, 0);
            for n_sensor = 1:N_sensor
                x_P{n_sensor,1} = ave_x_P(:,n_sensor);
                x_P{n_sensor,1} = N_sensor*x_P{n_sensor,1};
                x_P{n_sensor,1} = reshape(x_P{n_sensor,1},d,[]);
                bestEstimate = zeros(size(x_P{n_sensor,1},2),1);
                for n = 1:size(x_P{n_sensor,1},2)
                    [hp,~] = observation(x_P{n_sensor,1}(:,n),cd(:,n_sensor),0);
                    bestEstimate(n) = mvnpdf(hp', h_global{n_sensor,1}', R);
                end
    %             X = particle{n_sensor,1}; y = P_w{n_sensor,1}';
    %             [~, ~, ~, bestEstimate, ~] = GP_optfree(X, y, x, x_P{n_sensor,1}');
    %             P_w{n_sensor,1} = bestEstimate';
    %             data_weights = [data_weights P_w{n_sensor,1}'];
                data_weights = [data_weights bestEstimate];
            end
            ave_P_w = gossip(A, data_weights, iters_w, 0);  
            sum_P_w = ave_P_w*N_sensor;
            for n_sensor = 1:N_sensor
                P_w{n_sensor,1} = sum_P_w(:,n_sensor);
                [ x_P{n_sensor,1}, P_w{n_sensor,1} ] = Resample( x_P{n_sensor,1}, P_w{n_sensor,1}, N );
            end
            
        elseif strcmp(distributed, 'c')
            % 2.c) distributed resampling
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
                for n = 1:M/N_sensor
                    u =  rand; 
                    if u < p_var(1)
                        gm = gmdistribution(particle_set, sigmaf, pk); % based on alphaf
                    else
                        gm = gmdistribution(particle{n_sensor,1}, sigma, p); % based on alpha
                    end
                    x_P{n_sensor,1}(:,n) = random(gm,1).';
                end 

            end
            %%%%% communication (particle exchange)
            for n_sensor = 1:N_sensor 
                for neighbor = 1:length(neighbors{1,n_sensor})
                    x_P{neighbors{1,n_sensor}(neighbor),1} = ...
                        [x_P{neighbors{1,n_sensor}(neighbor),1} x_P{n_sensor ,1}(:,1:M/N_sensor)];
                end
            end

            for n_sensor = 1:N_sensor
                 P_w{n_sensor,1} = zeros(1,size(x_P{n_sensor,1},2)-M/N_sensor);
                 for n = M/N_sensor+1:size(x_P{n_sensor,1},2)
                    [hp,~] = observation(x_P{n_sensor,1}(:,n),cd(:,n_sensor),0);
                    P_w{n_sensor,1}(n-M/N_sensor) = mvnpdf(hp', h_global{n_sensor,1}', R);
                 end
                 tmp = x_P{n_sensor,1}(:,1:M/N_sensor);
                 tmp = tmp(:, randperm(M/N_sensor,round(M/N_sensor/(length(neighbors{1,n_sensor})+1))));
                 x_P{n_sensor,1}(:,1:M/N_sensor) = [];
                 [ x_P{n_sensor,1}, ~ ] = Resample( x_P{n_sensor,1}, P_w{n_sensor,1}, N-size(tmp,2) );
                 x_P{n_sensor,1} = [ tmp x_P{n_sensor,1}];
                 for n = 1:size(x_P{n_sensor,1},2)
                    x_P{n_sensor,1}(:,n) = x_P{n_sensor,1}(:,n) + chol(Sigma).'*randn(6,1);
                 end      
                 P_w{n_sensor,1} = ones(1,N)/N;
            end

            for iter = 1:iters
                for n_sensor = 1:N_sensor
                    x_P{n_sensor,1} = x_P{n_sensor,1}(:, randperm(N, round(M/N_sensor)));
                end
                for n_sensor = 1:N_sensor  
                    for neighbor = 1:length(neighbors{1,n_sensor})
                        x_P{neighbors{1,n_sensor}(neighbor),1} = ...
                            [x_P{neighbors{1,n_sensor}(neighbor),1} x_P{n_sensor ,1}(:,1:round(M/N_sensor))];
                    end
                end
                for n_sensor = 1:N_sensor
                    P_w{n_sensor,1} = zeros(1,size(x_P{n_sensor,1},2)-round(M/N_sensor));
                    for n = M/N_sensor+1:size(x_P{n_sensor,1},2)
                        [hp,~] = observation(x_P{n_sensor,1}(:,n),cd(:,n_sensor),0);
                        P_w{n_sensor,1}(n-M/N_sensor) = mvnpdf(hp', h_global{n_sensor,1}', R);
                    end
                    tmp = x_P{n_sensor,1}(:,1:M/N_sensor);
                    tmp = tmp(:, randperm(M/N_sensor,round(M/N_sensor/(length(neighbors{1,n_sensor})+1))));
                    x_P{n_sensor,1}(:,1:M/N_sensor) = [];
                    [ x_P{n_sensor,1}, ~ ] = Resample( x_P{n_sensor,1}, P_w{n_sensor,1}, N-size(tmp,2) );
                    x_P{n_sensor,1} = [ tmp x_P{n_sensor,1}];
                    for n = 1:size(x_P{n_sensor,1},2)
                        x_P{n_sensor,1}(:,n) = x_P{n_sensor,1}(:,n) + chol(Sigma).'*randn(6,1);
                    end
                    P_w{n_sensor,1} = ones(1,N)/N;
                end
            end

    %         for n_sensor = 1:N_sensor
    %              P_w{n_sensor,1} = zeros(1,size(x_P{n_sensor,1},2));
    %              for n = 1:size(x_P{n_sensor,1},2)
    %                 [hp,~] = observation(x_P{n_sensor,1}(:,n),cd(:,n_sensor),0);
    %                 P_w{n_sensor,1}(n) = mvnpdf(hp', h_global{n_sensor,1}', R);
    %              end
    %              [ x_P{n_sensor,1}, P_w{n_sensor,1} ] = Resample( x_P{n_sensor,1}, P_w{n_sensor,1}, N );
    %         end
    % 
    %         for iter = 1:iters
    %             for n_sensor = 1:N_sensor
    %                 x_P{n_sensor,1} = x_P{n_sensor,1}(:, randperm(N, ceil(M/N_sensor)));
    %             end
    %             for n_sensor = 1:N_sensor  
    %                 for neighbor = 1:length(neighbors{1,n_sensor})
    %                     x_P{neighbors{1,n_sensor}(neighbor),1} = ...
    %                         [x_P{neighbors{1,n_sensor}(neighbor),1} x_P{n_sensor ,1}];
    %                 end
    %             end
    %             for n_sensor = 1:N_sensor
    %                 P_w{n_sensor,1} = zeros(1,size(x_P{n_sensor,1},2));
    %                 for n = 1:size(x_P{n_sensor,1},2)
    %                     [hp,~] = observation(x_P{n_sensor,1}(:,n),cd(:,n_sensor),0);
    %                     P_w{n_sensor,1}(n) = mvnpdf(hp', h_global{n_sensor,1}', R);
    %                 end
    %                 [ x_P{n_sensor,1}, P_w{n_sensor,1} ] = Resample( x_P{n_sensor,1}, P_w{n_sensor,1}, N );
    %             end
    %         end
        end
    end
    
    
    % calculate estimated trajectory
    % use the data at the first sensor as a reference
    for n_sensor = 1:N_sensor
        x_est = mean(x_P{n_sensor,1},2); x_est_out{n_sensor,1} = [x_est_out{n_sensor,1} x_est];
    end
    error = 0;
    for n_sensor = 1:N_sensor
        error = error + norm(mean(x_P{n_sensor,1},2)-x);
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
%     figure;
%     plot(rmse);
%     title('RMSE Error');
%     xlabel('Time');ylabel('RMSE');
end

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
