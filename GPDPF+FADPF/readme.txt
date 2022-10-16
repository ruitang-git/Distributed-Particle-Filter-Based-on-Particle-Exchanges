--------------particle-based distributed particle filter------------

core function:

GP-DPF: particle_filter_GP2

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

NOTE: here iters_p = iters_w = [],and distributed = 'c';

FA-DPF: particle_filter_GPFA

function rmse = particle_filter_GPFA(N, iters_FA, iters_com, protocol, plt)

% Input
% N: particle size
% iters_FA: 1(no use)
% iters_com: number of communication iterations
% protocol: 'gossip'
% plt: plot figure or not 

To run this functions, please run run_sim.m;
Besides, the performance of all the DPF algorithms as well as the comparison results are included in this file;