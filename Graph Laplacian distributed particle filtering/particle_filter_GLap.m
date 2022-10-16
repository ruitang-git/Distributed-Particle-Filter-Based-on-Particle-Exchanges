function rms = particle_filter_GLap(plt,mode,num_sensors,topo,g_iters,m,GLA,if_subset)

% -------------------- %
% clear;clc;
% plt = 1;
% mode = 'multiple';
% num_sensors = 20;
% topo = 'distributed';
% g_iters = 2000;
% m = 100;
% GLA = 1;
% 
% for f = 1:10
% --------------------- %

%
% Description:
% 'This is core function of the reproduction work. It implements the
% distribued particle filter in the tracking problem.'
%  this file reproduces the work of paper 
%  "Rabbat, M., Coates, M., & Blouin, S. (2016, August). Graph Laplacian 
%  distributed particle filtering. In 2016 24th European Signal Processing 
%  Conference (EUSIPCO) (pp. 1493-1497). IEEE."
%
% INPUTS:
%  num_sensors - number of sensors
%  mode - decide to implement multi-sensor or single-sensor tracking 
%   = ~ 'multiple';
%       'single'
%  topo - underlyling communication topology
%   = ~ 'centralised'
%       'distributed'
%  GLA - decide the use of graph laplacian approximation    
%   = ~ 1 (yes)
%       0 (no)
%  g_iters - number of gossip iterations
%  m - approxiamtion (graph signal)
%  plt - decide the plotting process
%  if_subset - 0(default)

% OUTPUTS:
% rms - tracking error
%

configfile;

% observation of the initial states(variables)
% all variables with _ini are values at initialization step

% ---former initialization (abandon)---%
% theta_obsr = theta_ini+sqrt(theta_var)*randn; % theta_1 ~obtained by oberservation
% r_ini = r_mean+sqrt(r_var)*randn;
% s_ini = s_mean+sqrt(s_var)*randn; 
% c_mean = theta_obsr + pi/2;
% c_ini = c_mean+sqrt(c_var)*randn; % course angle
% x_true = [r_mean*sin(theta_obsr);r_mean*cos(theta_obsr);s_mean*sin(c_mean);s_mean*cos(c_mean)]; % true initial state
% x = [r_ini*sin(theta_ini);r_ini*cos(theta_ini);s_ini*sin(c_ini);s_ini*cos(c_ini)]; % estimated initial state vector(estimated from first measurement)
%

state = trajectory(STI,t,sigma_u);
x_ori = state(:,1); % origin state
d = size(state,1); % dimension of state vector
time_span = size(state,2) - 1; % do i need to minus 1 ?

% particle initialization
x_P = zeros(d,N); % particle
R = sigma_i^2 * eye(d);
for i = 1:N % particle sampling in 6-dimension space ? ? ?
    x_P(:,i) = x_ori + chol(R).' * randn(d,1);
end
    
% single-sensor tracking scenario
if strcmp(mode,'single')
    % sensor location
    l = [90;90];
    
    % varibles initialisation
    z_out = [];  % measurement value
    x_est_out = x_ori; % the vector of particle filter estimates.
    rms = []; % RMS position error
    
    % filtering process
    % G = [T^2/2 0;0 T^2/2;T 0;0 T];

    for sti = 1:time_span
% ----- former setting(abandon) ------- %
%         omega_true = am/sqrt(x_true(3)^2+x_true(4)^2);
%         F_true = [1 0 sin(omega_true*T)/omega_true -(1-cos(omega_true*T))/omega_true;...
%              0 1 (1-cos(omega_true*T))/omega_true sin(omega_true*T)/omega_true;...
%              0 0 cos(omega_true*T) -sin(omega_true*T);...
%              0 0 sin(omega_true*T) cos(omega_true*T)];
%         x_true = F_true * ( x_true + x_o ) - x_o + G * (repmat(v_mean',1,1) + randn(1,2)*sqrt(v_R))';
%         z = atan(x_true(1)/x_true(2)) + sqrt(w_var) * randn;
% ------------------------------------- %
        x_true =  state(:, sti+1);
        [z, R] = observation(x_true,l,sigma_v,sigma_w,1);
        x_P_update = zeros(d,N);
        for i = 1:N
            % sample particles from prior probability p(x(k)|x(k-1))
            
% ----- former setting(abandon) ------- %
%             omega = am/sqrt(x_P(3,i)^2+x_P(4,i)^2); % turning rate
%             F = [1 0 sin(omega*T)/omega -(1-cos(omega*T))/omega;...
%                  0 1 (1-cos(omega*T))/omega sin(omega*T)/omega;...
%                  0 0 cos(omega*T) -sin(omega*T);...
%                  0 0 sin(omega*T) cos(omega*T)];
%             x_P_update(:,i) = F * ( x_P(:,i) + x_o ) - x_o + G * (repmat(v_mean',1,1) + randn(1,2)*sqrt(v_R))';
%             % calculate the output value of particles, based on which to calculate the corresponding weights
%             z_update(i) = atan(x_P_update(1,i)/x_P_update(2,i));
% ------------------------------------- %
            x_P_update(:,i) = transition(x_P(:,i), t, sigma_u);
            [z_update(:,i),~] = observation(x_P_update(:,i),l,sigma_v,sigma_w,0);
            % calculate the weight of each particle, here we assume measurement
            % noise follows Gaussian distribution,
            % hereby w = p(y|x) can be calculated as follows
            P_w(i) = mvnpdf(z_update(:,i)', z', R);
        end
        % normalization.
        P_w = P_w./sum(P_w);

        % Resampling
        % also try histc function
        for i = 1 : N
            x_P(:,i) = x_P_update(:,find(rand <= cumsum(P_w),1));   % particles with higher weights get more offsprings
        end                                                    
        
        % state estimation, after resampling, the weights all become 1/N
        x_est = mean(x_P,2);

        error = norm(x_true-x_est);

        % Save data in arrays for later plotting
        z_out = [z_out z]; % observation
        x_est_out = [x_est_out x_est]; % estimated state
        rms = [rms error];
        
        % plot particle disribution
        if plt == 1
            figure;
            plot(x_P(1,:),x_P(2,:),'.',x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10);
        end
    end
      
    % plot 
    if plt == 1
        sti = 1:1+time_span;
        figure;
        subplot(2,2,1)
        plot(sti, state(1,:), '.-b', sti, x_est_out(1,:), '-.r','linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('x-coordinate');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('x-coordinate');
        subplot(2,2,2)
        plot(sti, state(2,:), '.-b', sti, x_est_out(2,:), '-.r','linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('y-coordinate');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('y-coordinate');
        subplot(2,2,3)
        plot(sti, state(3,:), '.-b', sti, x_est_out(3,:), '-.r', 'linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('x-velocity');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('x-velocity');
        subplot(2,2,4)
        plot(sti, state(4,:), '.-b', sti, x_est_out(4,:), '-.r','linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('y-velocity');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('y-velocity');
        figure;
        p1 = plot(state(1,:),state(2,:),'--o');
        hold on;
        p2 = plot(x_est_out(1,:),x_est_out(2,:),'--or'); % excludes the intial state estimate
        hold off;
        legend([p1,p2],{'true','estimated'});
        title('Agent trajectory');
    end

% multi-sensors scenario
elseif strcmp(mode,'multiple')
    % build sensor map
    [A, cd] = sensor_map(num_sensors, 200, 200, 0);
%     [A, cd] = sensor(num_sensors, roi);
%     x_o = [x_o; zeros(2,num_sensors)];
    
%     figure;
%     plot(state(1,:),state(2,:),'--o');
%     hold on;
%     gplot(A, cd.',':o');
%     hold on;
%     plot(cd(1,:),cd(2,:),'o','MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F');
%     legend('target trajectory','communication links','agents');
%     legend('Location','northeastoutside')
%     xticks(0:20:200);
%     yticks(0:20:200);
%     hold off;
    % varibles initialization
    z_out = zeros(num_sensors,1);
    x_est_out = x_ori; % the vector of particle filter estimates.
    rms = []; % RMS position error
    
    % filtering process
%     G = [T^2/2 0;0 T^2/2;T 0;0 T];
    for sti = 1:time_span
        % ----- former setting(abandon) ------- %
%         omega_true = am/sqrt(x_true(3)^2+x_true(4)^2);
%         F_true = [1 0 sin(omega_true*T)/omega_true -(1-cos(omega_true*T))/omega_true;...
%              0 1 (1-cos(omega_true*T))/omega_true sin(omega_true*T)/omega_true;...
%              0 0 cos(omega_true*T) -sin(omega_true*T);...
%              0 0 sin(omega_true*T) cos(omega_true*T)];
%         F_true = [1,0,1,0;0,1,0,1;0,0,1,0;0,0,0,1];
%         x_true = F_true *  x_true + G * (repmat(v_mean',1,1) + randn(1,2)*sqrt(v_R))';
        % -------------------------------------- %
        x_true = state(:, sti+1);
        z = zeros(2,num_sensors); % 2 observations are made
        for num_sensor = 1:num_sensors
            l = cd(:,num_sensor);
            [z(:,num_sensor),R] = observation(x_true,l,sigma_v,sigma_w,1);
        end      
        x_P_update = zeros(d,N);
        
        % only a subset(5) of sensor work(knn)
        est_state = transition(x_est_out(:,end),t,sigma_u);
        Idx = knnsearch(cd',est_state(1:2)','k',5);
        % y = randsample(n,k) 返回从整数 1 到 n 中无放回随机均匀抽取的 k 个值。
        
        % centralised
        if strcmp(topo, 'centralised')
            P_w = zeros(N,1);
            for i = 1:N
                % sample particles from prior probability p(x(k)|x(k-1))
                        % ----- former setting(abandon) ------- %
%                 omega = am/sqrt(x_P(3,i)^2+x_P(4,i)^2); % turning rate
%                 F = [1 0 sin(omega*T)/omega -(1-cos(omega*T))/omega;...
%                      0 1 (1-cos(omega*T))/omega sin(omega*T)/omega;...
%                      0 0 cos(omega*T) -sin(omega*T);...
%                      0 0 sin(omega*T) cos(omega*T)];
%                 F = [1,0,1,0;0,1,0,1;0,0,1,0;0,0,0,1];
%                 x_P_update(:,i) = F * x_P(:,i) + G * (repmat(v_mean',1,1) + randn(1,2)*sqrt(v_R))';
                       % -------------------------------------- %
                x_P_update(:,i) = transition(x_P(:,i),t,sigma_u);
                % calculate the output value of particles, based on which to calculate the corresponding weights
                for num_sensor = 1:num_sensors
                    l = cd(:,num_sensor);
                    [z_update(:,num_sensor),~] = observation(x_P_update(:,i),l,sigma_v,sigma_w,0);
                end
                % calculate the weight of each particle, here we assume measurement
                % noise follows Gaussian distribution,
                % hereby w = p(y|x) can be calculated as follows
%                 for sensor_k = 1:5 % only randomly choose 5 sensors to work
%                     num_sensor = Idx(sensor_k);
%                     P_w(i) = P_w(i)+ log(mvnpdf(z_update(:,num_sensor)', z(:,num_sensor)', R)); 
%                 end
                for num_sensor = 1:num_sensors 
                    P_w(i) = P_w(i)+ log(mvnpdf(z_update(:,num_sensor)', z(:,num_sensor)', R)); 
                end
                P_w(i) = exp(P_w(i));
            end
            
        % distributed
        elseif strcmp(topo, 'distributed')
            P_w = zeros(N,num_sensors);
            for i = 1:N
                % sample particles from prior probability p(x(k)|x(k-1))
                % ----- former setting(abandon) ------- %
%                 omega = am/sqrt(x_P(3,i)^2+x_P(4,i)^2); % turning rate
%                 F = [1 0 sin(omega*T)/omega -(1-cos(omega*T))/omega;...
%                      0 1 (1-cos(omega*T))/omega sin(omega*T)/omega;...
%                      0 0 cos(omega*T) -sin(omega*T);...
%                      0 0 sin(omega*T) cos(omega*T)];
%                 x_P_update(:,i) = F * x_P(:,i) + G * (repmat(v_mean',1,1) + randn(1,2)*sqrt(v_R))';
                % ------------------------------------- %
                x_P_update(:,i) = transition(x_P(:,i),t,sigma_u);
                % calculate the output value of particles, based on which to calculate the corresponding weights
                for num_sensor = 1:num_sensors
                    l = cd(:,num_sensor);
                    [z_update(:,num_sensor),~] = observation(x_P_update(:,i),l,sigma_v,sigma_w,0);
                end
                % calculate the weight of each particle, here we assume measurement
                % noise follows Gaussian distribution,
                % hereby w = p(y|x) can be calculated as follows
                if if_subset == 1
                    for sensor_k = 1:5 % only randomly choose 5 sensors to work
                        % y = randsample(n,k) 返回从整数 1 到 n 中无放回随机均匀抽取的 k 个值。
                        num_sensor = Idx(sensor_k);
                        P_w(i,num_sensor) = log(mvnpdf(z_update(:,num_sensor)', z(:,num_sensor)', R)); 
                    end
                else
                    for num_sensor = 1:num_sensors 
                        P_w(i,num_sensor) = log(mvnpdf(z_update(:,num_sensor)', z(:,num_sensor)', R)); 
                    end
                end
            end
            
            % GLA
            if GLA == 1
                
                % GLA_error_analysis;
                P_w(P_w < -30) = -30;
                [P_w_hat,V] = graphLap(x_P_update, N, P_w, m, k, 0); % approximation 
                
                % ----- check the approximation error ------ %
%                 P_w_pre = P_w; % save the data before approx.
%                 P_w = V*P_w_hat;
%                 error_GLA = norm(P_w-P_w_pre);
%                 P_w_c = gossip(A, P_w_pre, g_iters, 0);
%                 P_w_c = num_sensors * P_w_c;
%                 P_w_m = gossip(A,P_w_hat,g_iters, 0);
%                 P_w_d = V*P_w_m;
%                 P_w_d = num_sensors * P_w_d;
%                 error_GLA_1 = norm(P_w_d-P_w_c);
%                 fprintf('%f\n',min(min(P_w_pre)));
%                 debug_here = 1;
                
                % gossip
                P_w_hat = gossip(A,P_w_hat,g_iters, 0);
                P_w = V*P_w_hat;
                P_w = num_sensors * P_w; % summation
            else
                P_w = gossip(A, P_w, g_iters, 0);
                P_w = num_sensors * P_w; % summation
            end
            P_w = exp(P_w);        
        end
        
        % normalization.
        P_w = P_w./sum(P_w);

        if length(P_w(isnan(P_w)))
            % P_w(isnan(P_w)) = 1/length(P_w(isnan(P_w)));
            P_w = 1/N*ones(N,1);
        end
        
        % Resampling
        % also try histc function
        for i = 1 : N
            x_P(:,i) = x_P_update(:,find(rand <= cumsum(P_w),1));   % particles with higher weights get more offsprings
        end
        
        % state estimation, after resampling, the weights all become 1/N
        x_est = mean(x_P,2);

        error = norm(x_true-x_est);

        % Save data in arrays for later plotting
        z_out = [z_out z.']; % observation
        x_est_out = [x_est_out x_est]; % estimated state
        rms = [rms error];
        
        % plot particle disribution
%         if plt == 1
%             figure;
%             plot_particles(num_sensors, sti, x_P, 1/N*ones(N,1));
%             hold on;
%             plot(state(1,:),state(2,:),'--o')
%             hold off;
%         end   
%         if plt == 1
%             figure;
%             plot(x_P_update(1,:),x_P_update(2,:),'.',x_P(1,:),x_P(2,:),'.');
%             hold on;
%             plot(state(1,:),state(2,:),'--o');
%             hold on;
%             gplot(A, cd.',':o');
%             hold off;
%         end
    end
    % plot 
    if plt == 1
        sti = 1:1+time_span;
        figure;
        subplot(3,2,1)
        plot(sti, state(1,:), '.-b', sti, x_est_out(1,:), '-.r','linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('x-coordinate');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('x-coordinate');
        subplot(3,2,2)
        plot(sti, state(2,:), '.-b', sti, x_est_out(2,:), '-.r','linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('y-coordinate');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('y-coordinate');
        subplot(3,2,3)
        plot(sti, state(3,:), '.-b', sti, x_est_out(3,:), '-.r', 'linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('x-velocity');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('x-velocity');
        subplot(3,2,4)
        plot(sti, state(4,:), '.-b', sti, x_est_out(4,:), '-.r','linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('y-velocity');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('y-velocity');
        subplot(3,2,5)
        plot(sti, state(5,:), '.-b', sti, x_est_out(5,:), '-.r', 'linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('x-acceleration');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('x-acceleration');
        subplot(3,2,6)
        plot(sti, state(6,:), '.-b', sti, x_est_out(6,:), '-.r','linewidth',1);
        set(gca,'FontSize',12); set(gcf,'Color','White');
        xlabel('time step'); ylabel('y-acceleration');
        legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
        title('y-acceleration');
        figure;
        p1 = plot(state(1,:),state(2,:),'--o');
        hold on;
        p2 = plot(x_est_out(1,:),x_est_out(2,:),'--or'); % excludes the intial state estimate
        hold on;
        gplot(A, cd.',':o');
        hold on;
        plot(cd(1,:),cd(2,:),'o','MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F');
        hold off;
%         grid on;
        legend([p1,p2],{'true','estimated'});
%         title('Agent trajectory');
    end
end
% fprintf('%d finished!\n',f);
% end




