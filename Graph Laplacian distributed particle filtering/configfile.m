%%% Configuration file
%%% Permits various adjustments to parameters.

%% FORMER setting (abandon)
% % hyperparameter setup
% k = 20;
% % m = 100;
% 
% % parameters setup
% T = 1; % sampling time
% N = 1000; % number of particles
% 
% % state initialization
% r_mean = 10;
% r_var = 0.1; % range estimation distribution ~ Gaussian
% theta_ini = pi/4; % true initial bearing angle
% % theta_var = (1.5*(2*pi/360))^2;
% theta_var = 0.001;
% % s_mean = 0.3;
% s_mean = 0.05;
% s_var = 0; % speed distribution ~ Gaussian
% c_var = 0;
% particle_pos_var = 0.1; % particle distribution variance
% 
% % noise(process noise and measurement noise)
% % v_sigma2 = (1.6*10^(-6))^2;
% sigma = 1; % particle initialization
% v_sigma2 = 0.00035;
% v_mean = [0;0];
% v_var = v_sigma2*eye(2); 
% % v_R = chol(v_var); % process noise variance matrix (2x2)
% v_R = v_var;
% w_var = theta_var; % measurement noise variance
% am = 2*10^(-2); % manoeuvre acceleration
% 
% % agent network range
% len = 10; 
% width = 10;
%%
% hyperparameter setup
k = 5;
% m = 100;

% parameters setup
t = 1; % sampling time interval
sigma_u = 0.5;
sigma_v = 1;
sigma_w = 1;
sigma_i = 1; % particle initialization
STI = 30; % state transition intervals
N = 2000; % particle size
% N_sensor = 5; % number of sensors
roi = 60; % parameter determining the connectivity



