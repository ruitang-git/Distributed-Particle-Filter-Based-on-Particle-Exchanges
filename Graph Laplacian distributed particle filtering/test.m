close all;
clc;clear;
% configfile;

% parameters setup
T = 10; % sampling time
N = 10; % number of particles
rs = 0; % decide whether implement resampling
time_span = 50;

% noise(process noise and measurement noise)
v_sigma2 = 0.00035;
v_var = v_sigma2*eye(2); 
v_R = v_var;
w_var = 0.001; % measurement noise variance - attenuation
% CT Model
% w_var = (1.5*2*pi/360)^2; % bearing-angle
%
A = 10;
am = 1*10^(-2); % inversely proportional to the radius of curvature


% observation of the initial states(variables)
% all variables with _ini are values at initialization step

x_true = [4;4;0.05;0.05]; % true initial state
C = diag([2,2,0.001,0.001]);


% particle initialization
x_P = zeros(4,N); % particle
for i = 1:N 
    x_P(:,i) = x_true + sqrt(C) * randn(4,1);
end


% varibles initialisation
z_out = [];  % measurement value
x_est_out = []; % the vector of particle filter estimates.
true_state = [];
rms = []; % RMS position error

% filtering process
G = [T^2/2 0;0 T^2/2;T 0;0 T];

% CV Model
F = ...
    [1 0 T 0;...
     0 1 0 T;...
     0 0 1 0;...
     0 0 0 1];

for t = 1:time_span
    
    x_true = F * x_true + G * (sqrt(v_R)*randn(2,1));
    % CT Model
%     omega_true = am/sqrt(x_P(3,i)^2+x_P(4,i)^2); % turning rate
%     F_true = [1 0 sin(omega_true*T)/omega_true -(1-cos(omega_true*T))/omega_true;...
%          0 1 (1-cos(omega_true*T))/omega_true sin(omega_true*T)/omega_true;...
%          0 0 cos(omega_true*T) -sin(omega_true*T);...
%          0 0 sin(omega_true*T) cos(omega_true*T)];
%      
%     x_true = F_true * x_true + G * (sqrt(v_R)*randn(2,1));
    %
    
%     z = A/norm([x_true(1),x_true(2)]) + sqrt(w_var) * randn; % signal attenuation
    z = atan2(x_true(2),x_true(1)) + sqrt(w_var) * randn;
    
    x_P_update = zeros(4,N);
    for i = 1:N
        
        % CT Model
%         omega = am/sqrt(x_P(3,i)^2+x_P(4,i)^2); % turning rate
%         F = [1 0 sin(omega*T)/omega -(1-cos(omega*T))/omega;...
%              0 1 (1-cos(omega*T))/omega sin(omega*T)/omega;...
%              0 0 cos(omega*T) -sin(omega*T);...
%              0 0 sin(omega*T) cos(omega*T)];
         %
        % sample particles from prior probability p(x(k)|x(k-1))
        x_P_update(:,i) = F * x_P(:,i) + G * (sqrt(v_R)*randn(2,1));
        % calculate the output value of particles, based on which to calculate the corresponding weights
%         z_update(i) = A/norm([x_P_update(1,i),x_P_update(2,i)]);
        z_update(i) = atan2(x_P_update(2,i),x_P_update(1,i));
        % calculate the weight of each particle, here we assume measurement
        % noise follows Gaussian distribution,
        % hereby w = p(y|x) can be calculated as follows
        P_w(i) = 1/sqrt(2*pi*w_var) * exp(-(z - z_update(i))^2/(2*w_var)); 
    end
    % normalization.
    P_w = P_w./sum(P_w);

    % Resampling
    % also try histc function
    for i = 1 : N
        if rs == 1 % implement resample
            x_P(:,i) = x_P_update(:,find(rand <= cumsum(P_w),1));   % particles with higher weights get more offsprings
        else % not implement resample
            x_P(:,i) = x_P_update(:,i);
        end
    end                                                    


    % state estimation, after resampling, the weights all become 1/N
    if rs == 1
        x_est = mean(x_P,2);
    else 
        x_est = x_P * P_w';
    end

    error = norm(x_true(1:2)-x_est(1:2))^2;


    % Save data in arrays for later plotting
    true_state = [true_state x_true]; % true state
    % x_out = [x_out x];  %the actual output vector for measurement values.
    z_out = [z_out z]; % observation
    x_est_out = [x_est_out x_est]; % estimated state
    rms = [rms error];

    % plot particle disribution
    if  mod(t,time_span/5) == 1
%     if true 
        figure;
        plot_particles(1, t, x_P, 1/N*ones(N,1));
        hold on;
        plot(true_state(1,:),true_state(2,:),'--o')
        hold off;
    end
end



t = 1:time_span;
figure;
subplot(2,2,1)
plot(t, true_state(1,:), '.-b', t, x_est_out(1,:), '-.r','linewidth',1);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time step'); ylabel('x-coordinate');
legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
title('x-coordinate');
subplot(2,2,2)
plot(t, true_state(2,:), '.-b', t, x_est_out(2,:), '-.r','linewidth',1);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time step'); ylabel('y-coordinate');
legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
title('y-coordinate');
subplot(2,2,3)
plot(t, true_state(3,:), '.-b', t, x_est_out(3,:), '-.r', 'linewidth',1);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time step'); ylabel('x-velocity');
legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
title('x-velocity');
subplot(2,2,4)
plot(t, true_state(4,:), '.-b', t, x_est_out(4,:), '-.r','linewidth',1);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time step'); ylabel('y-velocity');
legend('true state', 'particle filter estimate','FontSize',8,'Location','northwest');
title('y-velocity');
figure;
plot(true_state(1,:),true_state(2,:),'-b')
hold on;
plot(x_est_out(1,:),x_est_out(2,:),'-r')
hold off;
legend('true trajectory','estimated trajectory');
title('Agent trajectory(CV Model)');


