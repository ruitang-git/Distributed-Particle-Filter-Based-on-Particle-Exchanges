%% configuration file
%
t = 1; % sampling time interval
sigma_u = 0.5;
sigma_v = 1;
sigma_w = 1;
sigma_i = 1; % particle initialization
STI = 30; % state transition intervals
N = 2000; % particle size
N_sensor = 9; % number of sensors
roi = 120; % parameter determining the connectivity
consensus_ite = 3;
f_tuning_ite = 3;
N_fusion = 5000;
N_recovery = 5000;
C = 1; % GMM components
maxiter = 100; % GMM iterations
reg = 1e-5;